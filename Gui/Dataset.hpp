#pragma once

#include <FL/Fl.H>
#include <Index.hpp>
#include <SimpleKmer.hpp>
#include <CsvIO.hpp>
#include <Centroid.hpp>
#include <FastaSequence.hpp>
#include <kNearestNeighbours.hpp>
#include <SparseSignature.hpp>
#include <SimpleKmer.hpp>

template<typename T>
class Dataset : public ICsvReader, public ICsvWriter {
protected:
	vector<T *> items;
	unordered_map<string, T *> index;
	const string tagName;
	function<T *()> constructor;
	function<void( T * t )> destructor;
	function<void( CsvReader & r, T & t )> readElement;
	function<void( CsvWriter & w, const T & t )> writeElement;
	function<const string &(const T & t)> getId;
	unordered_map<string, string> metadata;

public:
	Dataset() {}

	void Init(
		const string & tagName,
		decltype(constructor) constructor,
		decltype(destructor) destructor,
		decltype(readElement) readElement,
		decltype(writeElement) writeElement,
		decltype(getId) getId
	) {
		this->tagName = tagName;
		this->constructor = constructor;
		this->destructor = destructor;
		this->readElement = readElement;
		this->writeElement = writeElement;
		this->getId = getId;
	}

	virtual ~Dataset() {
		DestroyItems();
	}

	void DestroyItems() {
		for ( auto i : items ) destructor( i );
		items.clear();
	}

	T * Add( const T & element ) {
		auto &id = getId( element );
		auto &item = index[id];

		if ( item ) {
			throw Exception( "Duplicate Id in dataset!", FileAndLine );
		}

		auto t = constructor();
		*t = element;
		items.push_back( t );
		item = t;

		return t;
	}

	virtual void ParseMetadata() = 0;

	void Read( CsvReader & r ) {
		string id, name, value;
		size_t N;
		r >> id >> N;

		if ( id != tagName ) {
			throw Exception( "Expected tag " + tagName + " not found.", FileAndLine );
		}

		Clear();
		items.resize( N );

		while ( !r.IsEOF() ) {
			r >> id;

			if ( id != "meta" ) {
				r.Unget( id );
				break;
			}

			r >> name >> value;
			metadata[name] = value;
		}

		for ( size_t i = 0; (!r.IsEOF()) && (i < N); i++ ) {
			auto item = constructor();
			items[i] = item;
			readElement( r, *item );
			index[getId( *item )] = item;
		}
	}

	void Write( CsvWriter &w ) const {
		w << tagName << Size() << '\n';

		for ( auto & meta: metadata ) {
			w << "meta" << meta.first << meta.second << '\n';
		}

		for ( auto item : items ) {
			writeElement( w, *item);
		}
	}

	size_t Size() const { return items.size(); }

	T * operator[]( size_t i ) const { return items[i]; }

	T * operator[]( const string & id ) const {
		auto iter = index.find( id );
		return iter == index.end() ? nullptr : iter->second;
	}

	auto begin() const { return items.begin(); }

	auto end() const { return items.end(); }

	auto metaBegin() const { return metadata.begin(); }

	auto metaEnd() const { return metadata.end; }

	auto metaFind( const string & s ) const { return metadata.find( s ); }

	string & Metadata( const string & s ) {
		return metadata[s];
	}

	void Clear () {
		metadata.clear();
		index.clear();
		DestroyItems();
	}

	T * Back() const { return items.back(); }

	T * Front() const { return items.front(); }
};

class KmerDataset: public Dataset<SimpleKmer> {
private:
	const Dataset<FastaSequence> *sequences;
	uint kmerLength = 0;


	using T = SimpleKmer;
	static T * constructor() { return new SimpleKmer(); }

	static void destructor( T * t ) { delete t; }

	void ParseMetadata() {
		kmerLength = KmerLength();
	}

	void readElement( CsvReader & r, T & t ) {
		string id;
		uint position;
		r >> id >> position;
		auto seq = (*sequences)[id];
		SimpleKmer kmer_(seq, position, kmerLength);
		auto kmer = Add(s);

		while ( !(r.IsEOF() || r.IsEOL()) ) {
			r >> id >> position;
			kmer->Add( (*sequences)[id], position);
		}
	}

	static void writeElement( CsvWriter & w, const T & t ) {
		for ( auto & instance : t.Instances() ) {
			w << instance.sequence->IdStr() << instance.kmerPosition;
		}

		w.Ln();
	}

	const string &getId( const T & t ) {
		return t.IdStr();
	}
public:
	uint KmerLength() {
		return Uint::Parse(metadata[NAMEOF(KmerLength)]);
	}

	void SetKmerLength( uint value ) {
		kmerLength = value;
		metadata[NAMEOF( KmerLength )] = Util::ToString(value);
	}
};

class SequenceDataset : public Dataset<FastaSequence> {
private:

public:
};

class Vocabulary : public Dataset<Centroid> {
public:
};

class SignatureDataset : public Dataset<SparseSignature> {
	const SequenceDataset * sequences;
	const SimilarityMatrix * matrix;
	const vector<Centroid> * vocab;
	uint kmerLength;
	Distance thresholdKmerDistance;

public:
	using T = SparseSignature;

	SignatureDataset() {}

	void Init(
		const SequenceDataset * sequences,
		const SimilarityMatrix * matrix,
		const vector<Centroid> * vocab,
		uint kmerLength,
		Distance thresholdKmerDistance
	) {
		Dataset::Init( "signatures", constructor, destructor, [this]( CsvReader & r, T & t ) { readElement( r, t ); }, writeElement, getId );
		this->sequences = sequences;
		this->kmerLength = kmerLength;
		this->matrix = matrix;
		this->vocab = vocab;
		this->thresholdKmerDistance = thresholdKmerDistance;
	}

	static T *constructor() { return new SparseSignature(); }

	static void destructor( T * t ) { delete t; }

	void readElement( CsvReader & r, T & t ) {
		string id;
		r >> id;
		t.SetSequence( (*sequences)[id] );
		SparseSet &set = t;
		set.Read( r );
	}

	static void writeElement( CsvWriter & w, const T & t ) {
		throw NotImplementedException( FileAndLine );
	}

	static const string &getId( const T & t ) {
		static string empty;
		return t.Sequence() ? t.Sequence()->IdStr() : empty;
	}

};

