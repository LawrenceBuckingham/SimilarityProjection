#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <cstdio>
#include <cctype>

#include "Types.hpp"
#include "CharMap.hpp"
#include "PrecisionRecallRecord.hpp"
#include "Selector.hpp"
#include "Util.hpp"
#include "TrecEvalRecord.hpp"
#include "EncodedKmer.hpp"
#include "Alphabet.hpp"
#include <cstdint>
#include "CsvIO.hpp"
#include "Registry.hpp"
#include "SparseSet.hpp"
#include "Histogram.hpp"

namespace QutBio {

	/**
	 * Sequence data (currently rather coupled to the FASTA interchange format).
	 */
	class FastaSequence : public virtual ICsvWriter {
	protected:
		///	Normalised character representation of sequence. Normalisation is: 
		//	case conversion, alphabet encoding, and alphabet decoding.
		string charData;

		///  Sequence data stored as list of bytes.
		vector<Symbol> sequence;

		/// Bytes, recoded as a list of (overlapping) digrams. Thus, 
		//	digrams[i] = Horner(&bytes[i],alphabet->Size(),2);
		vector<Digram> digrams;

		///	Sequence metadata, parsed by splitting the definition line on '|' pipe
		/// symbols.
		vector<string> metadata;

		/// Zero-origin index within the metadata vector where the sequence ID
		///	may be found.
		size_t idIndex = 0;

		/// Should the char data be converted to lowercase when parsed. Default is
		///	true, to enable case-insensitive search and comparison. 
		bool forceLowerCase = true;

		///	The address of an Alphabet object which may be used to encode the sequence
		///	char data.
		Alphabet * alphabet = Alphabets::ByName( "DEFAULT" );

		/// Unique Id number, derived from metadata[idIndex];
		size_t idNumber = 0;

		/// Zero-origin (but optional) index within metadata vector of an entry containing class labels.
		/// Where multiple "classes" are present they are grouped within a single metadata field, separated by commas or semicolons.
		int classIndex = -1;

		/// Feature vector of canonical class reference numbers. These are obtained by splitting
		/// metadata[classIndex], then consulting a private static Registry to get canonical numbers.
		SparseSet classes;

	public:
		/**
		 * Construct a new FastaSequence from a string of sequence data.
		 *
		 * @param defLine The definition line which is presumed to contain metadata. Metadata
		 * 		is assumed to be organised into fields separated by '|' (pipe)
		 * 		symbols.
		 *
		 * @param sequence The sequence char data. This may be spread across multiple lines,
		 * 		and include spaces or '-' (minus) symbols. The minus symbols are assumed to be
		 * 		insertions created by alignment programs, and therefore ignored, as is any
		 * 		character ch for which isspace(ch) is true.
		 *
		 * @param idIndex The zero-origin index of the sequence ID in the metadata vector.
		 *
		 * @param alphabet The address of an alphabet which defines the set of possible symbols
		 * 		that may appear in the character data.
		 *
		 * @param forceLowerCase  Should the char data be converted to lowercase when parsed. Default is
		 *		true, to enable case-insensitive search and comparison.
		*/
		FastaSequence(
			const string &defLine,
			const string &sequence,
			int idIndex,
			Alphabet * alphabet,
			bool forceLowerCase = true
		) : idIndex( idIndex ), forceLowerCase( forceLowerCase ), alphabet( alphabet ) {
			SetSequence( sequence );
			SetDefLine( defLine );
		}

		/**
		 * Construct a new FastaSequence from a byte vector of sequence data.
		 *
		 * @param defLine The definition line which is presumed to contain metadata. Metadata
		 * 		is assumed to be organised into fields separated by '|' (pipe)
		 * 		symbols.
		 *
		 * @param sequence The sequence char data. This must not cantain anything but "real"
		 * 		sequence data. This is equivalent to the char data encoded by the alphabet.
		 *
		 * @param idIndex The zero-origin index of the sequence ID in the metadata vector.
		 *
		 * @param alphabet The address of an alphabet which defines the set of possible symbols
		 * 		that may appear in the character data.
		 *
		 * @param forceLowerCase  Should the char data be converted to lowercase when parsed. Default is
		 *		true, to enable case-insensitive search and comparison.
		*/
		FastaSequence(
			const string &defLine,
			const vector<Symbol> &sequence,
			int idIndex,
			Alphabet * alphabet,
			bool forceLowerCase = true
		) : sequence( sequence ), idIndex( idIndex ), forceLowerCase( forceLowerCase ), alphabet( alphabet ) {
			SetDefLine( defLine );
			charData = alphabet->Decode( sequence );
			alphabet->Encode( sequence.data(), sequence.size(), 2, digrams );
		}

		FastaSequence( const FastaSequence & other ) = delete;

		FastaSequence & operator=( const FastaSequence & other ) = delete;

		/**
		 *	Move constructor.
		 *	@param other A temporary object which is source for the move.
		 */
		FastaSequence( FastaSequence && other ) :
			charData( std::move( other.charData ) ),
			sequence( std::move( other.sequence ) ),
			metadata( std::move( other.metadata ) ),
			idIndex( other.idIndex ),
			idNumber( other.idNumber ) {
			this->digrams = std::move( other.digrams );
		}

		/**
		 *	Move assignment operator.
		 *	@param other A temporary object which is source for the move.
		 */
		FastaSequence & operator=( FastaSequence && other ) {
			sequence = std::move( other.sequence );
			metadata = std::move( other.metadata );
			charData = std::move( other.charData );
			idIndex = other.idIndex;
			idNumber = other.idNumber;
			digrams = other.digrams;
			return *this;
		}

		/**
		 * FastaSequence destructor.
		 */
		virtual ~FastaSequence() {}

		/**
		 * Get the ID of the sequence.
		 * @returns The ID of the sequence.
		 * @throws Exception is thrown if the idIndex equals or exceeds the
		 * 	size of the metadata vector.
		 */
		const string & IdStr() const {
			if ( idIndex >= metadata.size() ) {
				throw Exception( "Index Out Of Bounds: idIndex", FileAndLine );
			}

			return metadata[idIndex];
		}

		/**
		 *	Get the unique id number of this sequence.
		 */
		size_t Id() const {
			return idNumber;
		}

		/**
		 * Gets a readonly reference to the stored sequence byte data.
		 * @returns Reference to the byte sequence.
		 */
		const vector<Symbol> & Sequence() const {
			return sequence;
		}

		/**
		 * Ensures that the length of the sequence is at least `minLength`, padding with the 
		 * supplied default symbol if not.
		 */
		void EnsureLengthAtLeast( size_t minLength, Symbol defaultSymbol ) {
			string copy = charData;

			while ( copy.size() < minLength ) copy.push_back( alphabet->Decode(defaultSymbol) );

			SetSequence( copy );
		}

		/**
		 *	Gets a readonly reference to the stored sequence, encoded as a list
		 *	of overlapping digrams.
		 */
		const vector<Digram> & Digrams() const {
			return digrams;
		}

		/**
		 * Gets a copy of the char data, obtained by decoding the byte
		 * data with the alphabet.
		 * @returns A string containing a copy of the decoded byte data.
		 */
		const string & CharData() const {
			return charData;
		}

		/**
		 * Sets the byte data by encoding the supplied string character by
		 * character with the alphabet.
		 *
		 * @param value a string containing raw char data to encode. If
		 * 		forceLowerCase is true then the char data is converted to
		 * 		lower case prior to encoding. Symbol ch in value is discarded if
		 * 		ch = '-' or isspace(ch) is true.
		 */
		void SetSequence( const string &value ) {
			sequence.clear();

			for ( auto ch : value ) {
				if ( ch == '-' || isspace( ch ) )
					continue;

				if ( forceLowerCase ) ch = tolower( ch );

				auto byteVal = alphabet->Symbols().find(ch) == std::string::npos ? 
					alphabet->DefaultSymbol() :
					alphabet->Encode( ch );

				sequence.push_back( byteVal );
			}

			charData = alphabet->Decode( sequence );
			alphabet->Encode( sequence.data(), sequence.size(), 2, digrams );
		}

		/**
		 * Get a copy of the definition line, obtained by joining the metadata
		 * fields with '|' (pipe) separators.
		 *
		 * @returns A string obtained by concatenating the metadata fields with
		 * 		pipe symbols.
		 */
		string DefLine() const {
			return String::Join( metadata, "|" );
		}

		/**
		 * Gets a metadata value.
		 *
		 * @returns An element of the sequence metadata vector.
		 */
		const string & Metadata( int index ) const {
			static const string empty = "";
			return (0 <= index && (size_t) index < metadata.size()) ? metadata[index] : empty;
		}

		/**
		 *	Get the number of entries in the metadata vector.
		 */
		size_t MetaCount() const {
			return metadata.size();
		}

		/**
		 * Replaces the contents of the metadata vector by splitting the
		 * supplied string into fields separated by '|' (pipe) symbols.
		 *
		 * @param defLine The new definition line that will be applied to the
		 * 		sequence.
		 */

		void SetDefLine( const string &defLine ) {
			metadata = String::Split( defLine, '|' );

			if ( metadata.size() > 0 && metadata[0][0] == '>' ) {
				metadata[0].erase( 0, 1 );
			}

			idNumber = Register( IdStr() );
		}

		/**
		 * Gets the address of a dummy "zero" sequence which has no metadata and
		 * no char data.
		 *
		 * @returns The address of a static FastaSequence object.
		 */
		static FastaSequence * Zero() {
			static FastaSequence zero( "", "", 0, Alphabets::DEFAULT() );
			return &zero;
		}

		/// <summary>
		/// Parses a FASTA file and extracts a set of sequences. The fields of the
		/// definition line are assumed to be separated by '|' symbols, with an identifier
		/// in one field and the class label in another.
		/// <para>
		///		The id number may be a GI number, or something else such as a binding site ID or a
		///		application specific id number. This may be any string except one containing a pipe '|'
		///		symbol.
		/// </para>
		/// <para>
		///		The class label is an arbitrary string, which once again, cannot contain embedded pipe
		///		'|' symbols.
		/// </para>
		/// </summary>
		/// <param name="sequenceFile">The name of the file containing the sequences to be loaded.</param>
		/// <param name="idIndex">The index of the pipe-separated field (counting from 0) that contains the id number.</param>
		/// <param name="classIndex">The index of the pipe-separated field (counting from 0) that contains the class label.</param>
		/// <param name="sequences">A list onto which the sequences from the file will be appended.</returns>
		static void Read(
			istream &reader,
			int idIndex,
			Alphabet * alphabet,
			vector<FastaSequence *> & sequences
		) {
			string currentLine;
			string currentDefLine;
			string currentSequence;

			auto update = [&]() {
				if ( currentSequence.size() > 0 ) {
					sequences.push_back( new FastaSequence( currentDefLine, currentSequence, idIndex, alphabet ) );
				}
			};

			auto reset = [&]() {
				currentSequence.clear();
				currentDefLine = currentLine.substr( 1 );
			};

			while ( !(reader.fail() || reader.bad() || reader.eof()) ) {
				getline( reader, currentLine );
				String::TrimInPlace( currentLine );

				if ( currentLine[0] == '>' ) {
					update();
					reset();
				}
				else {
					currentSequence += String::Trim( currentLine );
				}
			}

			if ( !reader.eof() ) {
				throw Exception( "Error reading from stream.", __FILE__, __LINE__ );
			}

			update();

			// Console.Error( String.Format( "{0} definition lines read.", sequences.Count ) );
		}

		/// <summary>
		/// Parses a FASTA file and extracts a set of sequences. The fields of the
		/// definition line are assumed to be separated by '|' symbols, with an identifier
		/// in one field and the class label in another.
		/// <para>
		///		The id number may be a GI number, or something else such as a binding site ID or a
		///		application specific id number. This may be any string except one containing a pipe '|'
		///		symbol.
		/// </para>
		/// <para>
		///		The class label is an arbitrary string, which once again, cannot contain embedded pipe
		///		'|' symbols.
		/// </para>
		/// </summary>
		/// <param name="sequenceFile">The name of the file containing the sequences to be loaded.</param>
		/// <param name="idIndex">The index of the pipe-separated field (counting from 0) that contains the id number.</param>
		/// <param name="classIndex">The index of the pipe-separated field (counting from 0) that contains the class label.</param>
		/// <param name="sequences">A list onto which the sequences from the file will be appended.</returns>
		static vector<FastaSequence *> Read(
			istream &reader,
			int idIndex,
			Alphabet * alphabet
		) {
			vector<FastaSequence *> db;
			Read( reader, idIndex, alphabet, db );
			return db;
		}

		/// <summary>
		/// Parses a FASTA file and extracts a set of sequences. The fields of the
		/// definition line are assumed to be separated by '|' symbols, with an identifier
		/// in one field and the class label in another.
		/// <para>
		///		The id number may be a GI number, or something else such as a binding site ID or a
		///		application specific id number. This may be any string except one containing a pipe '|'
		///		symbol.
		/// </para>
		/// <para>
		///		The class label is an arbitrary string, which once again, cannot contain embedded pipe
		///		'|' symbols.
		/// </para>
		/// </summary>
		/// <param name="sequenceFile">The name of the file containing the sequences to be loaded.</param>
		/// <param name="idIndex">The index of the pipe-separated field (counting from 0) that contains the id number.</param>
		/// <param name="classIndex">The index of the pipe-separated field (counting from 0) that contains the class label.</param>
		/// <returns>A dictionary containing the sequences from the file, indexed by their id numbers.</returns>

		static vector<FastaSequence *> Read(
			const string &fileName,
			int idIndex,
			Alphabet * alphabet
		) {
			ifstream reader( fileName );
			return Read( reader, idIndex, alphabet );
		}

		static void Read(
			vector<FastaSequence *> & result,
			const string &fileName,
			int idIndex,
			Alphabet * alphabet
		) {
			ifstream reader( fileName );
			Read( reader, idIndex, alphabet, result );
		}

		static void Read(
			vector<FastaSequence *> & db,
			CsvReader & reader,
			int idIndex,
			Alphabet * alphabet,
			bool forceLowerCase = true
		) {
			db.clear();
			FastaSequence *seq = 0;

			while ( 0 != (seq = Read( reader, idIndex, alphabet, forceLowerCase )) ) {
				db.push_back( seq );
			}
		}

		/// <summary> Serialises this FastaSequence to text.
		/// </summary>
		/// <returns></returns>

		string ToString() {
			ostringstream out;
			out << (this);
			return out.str();
		}

		/**
		 *	<summary>
		 *		Gets a normalised histogram with the counts of each symbol.
		 *	</summary>
		 *	<param name="db">
		 *		Reference to a pointer-list of FastaSequence objects in which symbol occurrences are to be calculated.
		 *	</param>
		 *	<param name="histogram">
		 *		Reference to a histogram that will be overwritten with the symbol probabilities observed in the database.
		 *	</param>
		 */

		static Histogram<Symbol> GetSymbolHistogram( const vector<FastaSequence *> & db ) {
			Histogram<Symbol> histogram;

			for ( auto seq : db ) {
				histogram.AddRange( seq->sequence );
			}

			histogram.Normalise();
			return histogram;
		}

		/**********************
		 *	<summary>
		 *		Uses the supplied selector to choose zero or more kmers from this sequence and pass them back for processing.
		 *	</summary>
		 *	<param name="kmerLength">The size of the desired kmers.</param>
		 *	<param name="selector">A reference to a Selector that determines which, if any, kmers are processed.</param>
		 *	<param name="process">A function that is invoked once per selected kmer to process the selection.</param>
		 **********************/

		void SelectKmers(
			size_t kmerLength,
			Selector &selector,
			function<void( FastaSequence *seq, size_t pos, size_t length )> process
		) {
			size_t length = sequence.size();

			if ( length < kmerLength )
				return;

			size_t n = length - kmerLength + 1;

			for ( size_t i = 0; i < n; i++ ) {
				if ( selector.SelectThis() ) {
					process( this, i, kmerLength );
				}
			}
		}

		/**
		 *	<summary>
		 *		Iterates over the kmers in this sequence and passes them back for processing.
		 *	</summary>
		 *	<param name="kmerLength">The size of the desired kmers.</param>
		 *	<param name="process">A function that is invoked once per selected kmer to process the selection.</param>
		 */

		void SelectKmers(
			size_t kmerLength,
			function<void( FastaSequence * seq, size_t pos, size_t length )> process
		) {
			size_t n = KmerCount( kmerLength );

			for ( size_t i = 0; i < n; i++ ) {
				process( this, i, kmerLength );
			}
		}

		/**
		*	<summary>
		*		Returns the total number of kmers in a dataset.
		*	</summary>
		*	<param name="db">A list of sequences.</param>
		*	<param name="kmerLength">The kmer length.</param>
		*/
		static size_t GetTotalKmerCount(
			const vector<FastaSequence *> & db,
			size_t kmerLength
		) {
			size_t totalKmerCount = 0;

			for ( auto seq : db ) {
				totalKmerCount += seq->KmerCount( kmerLength );
			}

			return totalKmerCount;
		}

		/**
		 *	Returns the number of complete kmers of length K in the sequence.
		 *	In particular, if there are fewer than K symbols in the sequence, it
		 *	returns 0.
		 *	@param K the word (k-mer) length.
		 *	@param Returns max(0, length - K + 1).
		 */
		size_t KmerCount( size_t K ) const {
			auto length = sequence.size();
			return length >= K ? length + 1 - K : 0;
		}

		/**
		 * Gets the length of the sequence.
		 * @returns The length of the sequence.
		 */
		size_t Length() const {
			return sequence.size();
		}

		/**
		 * Serialises the sequence to a stream in FASTA format, as two lines containing
		 * definition (concatenated metadata) and sequence respectively. Note that
		 * a linefeed is emitted after the sequence.
		 *
		 * @param out Reference to an output stream.
		 * @param seq The address of a sequence which is to be serialised.
		 * @returns The stream.
		 */

		friend ostream &operator<<( ostream &out, const FastaSequence * seq ) {
			out << '>' << seq->DefLine() << endl;

			for ( auto c : seq->sequence ) {
				out << seq->alphabet->Decode( c );
			}

			out << endl;
			return out;
		}

		/**
		 * Serialises the current sequence to a stream in FASTA format, as two lines containing
		 * definition (concatenated metadata) and sequence respectively. Note that
		 * a linefeed is emitted after the sequence.
		 *
		 * @param out Pointer to a FILE struct which encapsulates the stream.
		 */

		void fprint( FILE * out ) {
			fprintf( out, ">%s\n%s\n", DefLine().c_str(), CharData().c_str() );
		}

		/**
		 *	Serialise this FASTA sequence to a CSV-formatted document.
		 */
		void Write( CsvWriter &w ) const override {
			w.WriteFields( metadata ).Sep().Write( CharData() ).Ln();
		}

		static FastaSequence * Read( CsvReader & csvReader, int idIndex, Alphabet * alphabet, bool forceLowerCase = true ) {
			FastaSequence *seq = new FastaSequence( "", "", idIndex, alphabet, forceLowerCase );

			csvReader.ReadRecord( seq->metadata );

			if ( seq->metadata.size() < 2 ) {
				delete seq;
				return 0;
			}

			seq->SetSequence( seq->metadata.back() );
			seq->metadata.pop_back();
			return seq;
		}

		/**
		 *	<summary>
		 *		Pads the sequence out to the designated minimum length
		 *	</summary>
		 *	<param name="minLength">The required minimum length.</param>
		 *	<param name="padding">The symbol used to pad the sequence.</param>
		 */

		void Pad( size_t minLength, Symbol padding ) {
			while ( sequence.size() < minLength ) {
				sequence.push_back( padding );
			}
		}

		/**
		 *	Count how many occurrences there are of each distinct metadata entry (by index).
		 *
		 *	@param db A collection of FastaSequence addresses to be surveyed.
		 *	@param metaCount A list of mappings which will be used as bags to count distinct metadata values. PRevisous contents will be overwritten.
		 *	@pre seq in db ==> seq != null.
		 *	@post metaCount.size() = db[0].metadata.size()
				and (dom(metaCount[i]) = {seq:db --> seq->metadata[i]})
				and (s in dom(metaCount[i]) ==> metaCount[i][s] = #{seq: db | seq->metadata[i] = s})
				)$
		 */
		static void GetMetadataCounts( const vector<FastaSequence *> & db, vector<map<string, size_t>> & metaCount ) {
			metaCount.resize( db[0]->metadata.size() );

			for ( auto seq : db ) {
				auto &metadata = seq->metadata;

				for ( size_t i = 0; i < metadata.size(); i++ ) {
					metaCount[i][metadata[i]] ++;
				}
			}
		}

		/**
		 *	Get canonical internal numeric id number from id string.
		 *	@param idStr A sequence id tag.
		 *	@returns A distinct number which corresponds to the id tag.
		 */
		static size_t Register( const string & idStr ) {
			static Registry<string, std::hash<string>> register_;
			return register_( idStr );
		};

		/**
		 *	Get canonical internal numeric id number from id string.
		 *	@param classLabel A class id tag.
		 *	@returns A distinct number which corresponds to the class id tag.
		 */
		static size_t RegisterClass( const string & classLabel ) {
			static Registry<string, std::hash<string>> register_;
			return register_( classLabel );
		};

		const SparseSet & Classes() const { return classes; }

		int ClassIndex() const { return classIndex; }

		void SetClassIndex( int value ) {
			if ( ClassIndex() == value ) return;

			classIndex = value;
			classes.Clear();

			if ( value < 0 || value >= (int) metadata.size() ) return;

			string classLabel;
			istringstream istr( metadata[classIndex] );
			CsvReader r( istr, ';' );

			while ( !r.IsEOF() ) {
				r >> classLabel;
				classes.Add( RegisterClass( classLabel ) );
			}

			classes.Sort();
		}

		/**
		 *	Determine whether this sequence is related to another, in the sense
		 *	that their respective lists of class labels have non-empty intersection.
		 *	@param other The address of a sequence for comparison.
		 *	@returns true iff
		 */
		bool IsRelated( const FastaSequence * other ) const {
			return classes.Similarity( other->classes ) > 0;
		}

	private:
		int nameIndex = -1;

	public:
		/// Get the optional position of the sequence name in the metadata list.
		int NameIndex() const { return nameIndex; }

		/**
		 *	Set the  optional position of the sequence name in the metadata list.
		 *	@param nameIndex	The zero-origin position of the sequence name within
		 *						the metadata list. Indexes outside dom(metadata) are
		 *						treated as "no name".
		 */
		void SetNameIndex( int nameIndex ) { this->nameIndex = nameIndex; }

		/**
		 *	Get the name of the sequence.
		 *	@returns nameIndex in dom(metadata) ? metadata[nameIndex] : "".
		 */
		const string & Name() const {
			static string empty = "";
			return 0 <= nameIndex && nameIndex < (int) metadata.size() ? metadata[nameIndex] : empty;
		}
	};

	/**
	 * Symbolic name for "Pointer to FastaSequence" type.
	 */
	using pFastaSequence = FastaSequence * ;
} // namespace QutBio
