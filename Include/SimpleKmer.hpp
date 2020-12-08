#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <ios>
#include <unordered_map>
#include <unordered_set>

#include "FastaSequence.hpp"
#include "SubString.hpp"
#include "DistanceType.hpp"

namespace QutBio {

	using pSimpleKmer = class SimpleKmer *;

	class SimpleKmer :
		public virtual ICsvWriter {
	public:
		/*
		** <summary>The location of an instance of the Kmer.</summary>
		*/
		struct Instance :
			public virtual ICsvWriter {
			// The Id of the containing sequence.
			const FastaSequence * sequence;

			// The offset of the first character of the kmer from the start of the sequence.
			size_t kmerPosition;

			Instance( const FastaSequence * sequence, size_t kmerPosition ) : sequence( sequence ), kmerPosition( kmerPosition ) {}

			const Symbol * Bytes() const {
				return sequence->Sequence().data() + kmerPosition;
			}

			const Digram * Digrams() const {
				return sequence->Digrams().data() + kmerPosition;
			}

			friend ostream & operator<<( ostream & str, const Instance & instance ) {
				str << instance.sequence->IdStr() << ":" << instance.kmerPosition;
				return str;
			}

			static Instance & Zero() {
				static Instance zero( FastaSequence::Zero(), 0 );
				return zero;
			}

			friend bool operator==( const Instance & lhs, const Instance & rhs ) {
				return lhs.sequence == rhs.sequence && lhs.kmerPosition == rhs.kmerPosition;
			}

			friend bool operator!=( const Instance & lhs, const Instance & rhs ) {
				return lhs.sequence != rhs.sequence || lhs.kmerPosition != rhs.kmerPosition;
			}

			void Write( CsvWriter &w ) const override {
				w.Write( sequence->IdStr() ).Write( kmerPosition );
			}

			string Chars( size_t kmerLength ) const {
				return sequence->CharData().substr( kmerPosition, kmerLength );
			}

		};

	private:
		// The string representation of the kmer.
		Substring substring;

		// The locations of all instances of the current kmer.
		vector<Instance> instances;

		size_t serialNumber = 0;

	public:
		/**
		**	Summary:
		**		Construct a kmer belonging to a sequence.
		**
		**	Parameters:
		**		seq:          a sequence containing a prototypical instance of the kmer;
		**		kmerPosition: the offset of the kmer from the start of the sequence;
		**		kmerLength:   the length of the kmer;
		**		dist:         a distance value associated with the kmer, used when clustering
		**			          Arguably this should not be present in the base class.
		*/
		SimpleKmer(
			FastaSequence * seq,
			size_t kmerPosition,
			size_t kmerLength
		) :
			substring( seq->Sequence().data(), kmerPosition, kmerLength )
			//
		{
			Add( seq, kmerPosition );
		}

		virtual ~SimpleKmer() {

		}

		/**
		 *	Gets the byte array associated with the first instance of the kmer.
		 *	@returns FirstInstance().Bytes(). If there is no first instance, then this will return the null pointer encapsulated by Instance::Zero.
		 */
		const Symbol * Bytes() const {
			return FirstInstance().Bytes();
		}

		/**
		**	Summary:
		**		Adds a new instance to the current kmer.
		**
		**	Parameters:
		**		seq:          a sequence containing a prototypical instance of the kmer;
		**		kmerPosition: the offset of the kmer from the start of the sequence;
		**		kmerLength:   the length of the kmer;
		**		dist:         a distance value associated with the kmer, used when clustering
		**			          Arguably this should not be present in the base class.
		*/
		void Add( FastaSequence * seq, size_t kmerPosition ) {
			instances.emplace_back( seq, kmerPosition );
		}

		template<typename CollectionType>
		void Add( const CollectionType & other ) {
			for ( auto & i : other ) {
				Add( i.sequence, i.kmerPosition );
			}
		}

		const Substring & Substr() const {
			return substring;
		}

		/// Returns a vector of byte containing a copy of the kmer.
		///	@returns A vector containing a copy of the bytes in this substring.

		vector<Symbol> Word() const {
			vector<Symbol> ret( substring.Chars(), substring.Chars() + substring.Length() );
			return ret;
		}

		///	Get the normalised character representation of the substring.
		///	@returns If instance list is empty, returns a string of length zero. Otherwise, returns a new string containing the normalised substring of length k starting at the k-mer index of the first instance. 

		string Chars() const {
			if ( instances.size() == 0 ) return "";

			auto i = instances.front();
			return i.Chars( this->Length() );
		}

		/// Gets a readonly reference to the instance list for this k-mer.

		const vector<Instance> & Instances() const {
			return instances;
		}

		friend ostream & operator << ( ostream & str, const SimpleKmer & kmer ) {
			for ( auto & instance : kmer.instances ) {
				str << instance << ";";
			}

			return str;
		}

		Kmer & operator=( const SimpleKmer & other ) = delete;

		friend bool operator==( const SimpleKmer & lhs, const SimpleKmer & rhs ) {
			return lhs.substring == rhs.substring;
		}

		friend bool operator!=( const SimpleKmer & lhs, const SimpleKmer & rhs ) {
			return lhs.substring != rhs.substring;
		}

		friend bool operator<( const SimpleKmer & lhs, const SimpleKmer & rhs ) {
			return lhs.substring < rhs.substring;
		}

		ostream & write( ostream & stream ) const {
			return stream << substring;
		}

		/**
		 *	<summary>
		 *		Returns the number of kmers required to tile the longest sequence in a dataset.
		 *	</summary>
		 *	<param name="db">A list of sequences.</param>
		 *	<param name="kmerLength">The kmer length.</param>
		 */
		static size_t GetMaxKmerCount( vector<FastaSequence *> & db, size_t kmerLength ) {
			size_t maxKmerCount = 0;

			for ( auto seq : db ) {
				size_t K = seq->KmerCount( kmerLength );

				if ( K > maxKmerCount ) {
					maxKmerCount = K;
				}
			}

			return maxKmerCount;
		}

		const FastaSequence * Sequence() const {
			return instances.size() > 0 ? instances[0].sequence : FastaSequence::Zero();
		}

		size_t KmerPosition() const { return instances.size() > 0 ? instances[0].kmerPosition : 0; }

		const Instance & FirstInstance() const { return instances.size() > 0 ? instances[0] : Instance::Zero(); }

		size_t Length() const { return substring.Length(); }

		size_t SerialNumber() const { return serialNumber; }

		SimpleKmer & SetSerialNumber( size_t value ) {
			serialNumber = value;
			return *this;
		}

		string Meta( int index ) {
			return instances.size() == 0 ? string( "" ) : instances[0].sequence->Metadata( index );
		}

		void Write( CsvWriter & writer ) const override {
			writer.Write( FirstInstance() );
		}

		struct Subsequence {
			FastaSequence * source;
			size_t start;
			size_t length;
		};

		/// Map Substring to SimpleKmer containing exact instances of the Substring. 
		class Index : public unordered_map<Substring, SimpleKmer *, Substring::Hash> {
			// TODO: this is not a valid use of inheritance. Convert the map to a member variable.

		public:
			using BaseType = unordered_map<Substring, SimpleKmer *, Substring::Hash>;

			Index() {}

			/*
				**	Summary:
				**		Constructs a KmerIndex.
				**
				**	Parameters:
				**		dataset:    a list of sequences;
				**		kmerLength: the word length for tiling;
				*/

			Index(
				const vector<FastaSequence *> & dataset,
				size_t kmerLength
			) {
				AddRange( dataset, kmerLength );
			}

			void AddRange( const std::vector<QutBio::FastaSequence *> & dataset, const size_t &kmerLength ) {
				size_t length = dataset.size();

				for ( size_t i = 0; i < length; i++ ) {
					auto seq = dataset[i];
					const Symbol *residues{ seq->Sequence().data() };
					const size_t seqLen{ seq->Sequence().size() };

					if ( seqLen < kmerLength ) continue;

					uint kmerCount = seq->KmerCount( kmerLength );

					for ( size_t kmerPos = 0; kmerPos < kmerCount; kmerPos++ ) {
						AddKmer( residues, kmerPos, kmerLength, seq );
					}
				}
			}

			Index(
				const vector<FastaSequence *> & dataset,
				const vector<size_t> & selection,
				size_t kmerLength
			) {
				AddRange( dataset, selection, kmerLength );
			}

			void AddRange(
				const std::vector<QutBio::FastaSequence *> & dataset,
				const std::vector<size_t> & selection,
				size_t kmerLength
			) {
				size_t length = selection.size();

				for ( size_t i = 0; i < length; i++ ) {
					auto seq = dataset[selection[i]];
					Add( seq, kmerLength );
				}
			}

			void Add(
				FastaSequence *seq,
				size_t kmerLength
			) {
				const Symbol *residues{ seq->Sequence().data() };
				const size_t seqLen{ seq->Sequence().size() };
				uint kmerCount = seq->KmerCount( kmerLength );

				for ( size_t kmerPos = 0; kmerPos < kmerCount; kmerPos++ ) {
					AddKmer( residues, kmerPos, kmerLength, seq );
				}
			}

			Index( const Index && other ) : BaseType( other ) {}

			Index & operator= ( const Index && other ) {
				BaseType::operator=( other );
				return *this;
			}

			void AddKmer(
				const Symbol * residues,
				const size_t kmerPos,
				const size_t kmerLength,
				FastaSequence * seq
			) {
				Substring key( residues, kmerPos, kmerLength );
				auto item = this->find( key );

				if ( item == this->end() ) {
					auto value = new SimpleKmer( seq, kmerPos, kmerLength );
					std::pair<Substring, SimpleKmer *> p( key, value );
					BaseType::insert( p );
				}
				else {
					item->second->Add( seq, kmerPos );
				}
			}

			Index(
				vector<Subsequence> & subSequences,
				size_t kmerLength
			) {
				for ( auto & subSequence : subSequences ) {
					auto seq = subSequence.source;
					const Symbol *residues{ seq->Sequence().data() };
					uint kmerCount = seq->KmerCount( kmerLength );

					for ( size_t kmerPos = subSequence.start;
						kmerPos < kmerCount && kmerPos + kmerLength <= subSequence.start + subSequence.length;
						kmerPos++
						) {
						AddKmer( residues, kmerPos, kmerLength, seq );
					}
				}
			}

			Index(
				const Subsequence * subSequences,
				size_t length,
				size_t kmerLength
				//
			) {
				for ( size_t i = 0; i < length; i++ ) {
					auto subSequence = subSequences[i];
					auto seq = subSequence.source;
					auto residues{ seq->Sequence().data() };
					auto kmerCount = seq->KmerCount( kmerLength );

					for ( auto kmerPos = subSequence.start; kmerPos < kmerCount && kmerPos + kmerLength <= subSequence.start + subSequence.length; kmerPos++ ) {
						AddKmer( residues, kmerPos, kmerLength, seq );
					}
				}
			}

			/**
				**	Summary:
				**		Destructor.
				*/
			virtual ~Index() {
				for ( auto &p : *this ) {
					delete p.second;
				}
			}

			/**
				**	Summary:
				**		Updates an internal cache of kmer pointers and then returns it by reference.
				*/
			void GetKmers( vector<SimpleKmer *> & allKmers ) const {
				size_t id = 0;

				allKmers.clear();

				for ( auto &pair : *this ) {
					auto kmer = pair.second;
					kmer->SetSerialNumber( id++ );
					allKmers.push_back( kmer );
				}
			}

			void GetInstances( vector<Instance> &instances ) const {
				instances.clear();

				for ( auto &pair : *this ) {
					auto kmer = pair.second;
					for ( auto instance : kmer->instances ) {
						instances.push_back( instance );
					}
				}
			}

			size_t size() {
				return BaseType::size();
			}
		};

		static const SimpleKmer &Zero() {
			static SimpleKmer zero( FastaSequence::Zero(), 0, 0 );
			return zero;
		}
	};
}

