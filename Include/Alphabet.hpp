#pragma once

#include <string>
#include <cstdint>
#include <cstddef>
#include <map>

#include "Assert.hpp"
#include "EncodedKmer.hpp"
#include "EnumBase.hpp"
#include "Array.hpp"
#include "BitSet.hpp"

using namespace QutBio;

namespace QutBio {
	typedef class Alphabet *pAlphabet;

	class Alphabet {

	protected:
		string name;
		string symbols;
		vector<byte> inverse;
		Symbol defaultSymbol;

		/// <summary>Define the default symbol for the alphabet.</summary>
		void DefaultSymbol( Symbol symbol ) {
			defaultSymbol = symbol;
		}

	public:
		Alphabet( string name, string symbols )
			: name( name ), symbols( symbols ), inverse( 256 ), defaultSymbol( Symbol::From(0) ) {

			std::fill( inverse.begin(), inverse.end(), 0 );

			for ( size_t i = 0; i < symbols.size(); i++ ) {
				inverse[tolower( symbols[i] )] = (unsigned char) i;
				inverse[toupper( symbols[i] )] = (unsigned char) i;
			}
		}

		/**
		 *	<summary>
		 *	Computes the number of KmerWord values required to pack a sequence of bytes
		 *	taken charsPerWord at a time.
		 *	</summary>
		 *	<param name="length"></summary>
		 *	<param name="length"></summary>
		 *	<returns></returns>
		 */
		static size_t WordsRequiredToPack( size_t length, size_t charsPerWord ) {
			return (length + charsPerWord - 1) / charsPerWord;
		}

		/** <summary>
		 *	Packs a sequence (in my initial formulations, this was assumed to be a
		 *	single kmer) as a sequence of numeric values, inserting ${charsPerWord} items
		 *	in each word of the output array.
		 * </summary>
		 * <param name="s">An array of bytes (unsigned 8 bit integers) to be packed.</param>
		 * <param name="charsPerWord">
		 *	The number of encoded characters to pack into each word.
		 *	This is independent of any particular kmer tiling scheme. It has to do with the
		 *	maximum size of a pre-computed kmer similarity lookup table that can be created.
		 *	<para>
		 *		For my current protein experiments, I use charsPerWord == 2 or sometimes 3
		 *		for amino acid, with the BLOSUM alphabet, this yields tables of approx 1/4GB.
		 *	</para>
		 * </param>
		 * <returns></returns>
		 */

		virtual void Encode( const Symbol s[], size_t kmerLength, size_t charsPerWord, EncodedKmer code ) const {
			size_t words = Alphabet::WordsRequiredToPack( kmerLength, charsPerWord );
			size_t size = symbols.size();
			size_t wordIndex = 0;
			code[wordIndex] = 0;

			for ( size_t i = 0; i < kmerLength; i++ ) {
				int j = s[i].value;
				code[wordIndex] = code[wordIndex] * (int) size + j;

				if ( i > 0 && i % charsPerWord == (charsPerWord - 1) ) {
					wordIndex++;

					if ( wordIndex < words ) {
						code[wordIndex] = 0;
					}
				}
			}
		}

		/**
		 *	Alternative experimental coding which packs di-grams or tri-grams into
		 *	a list of some suitable unsigned integer data type.
		 *	@typeparam	T Probably uint16_t, or something like that.
		 *	@param[in]	s A sequence of bytes to be converted into a sequence of packed n-grams.
		 *	@param[in]	len The number of symbols in s to process.
		 *	@param[in]	charsPerWord The 'n' in n-gram.
		 *	@param[out]	code A vector such that code[offset]=Horner(s[offset],size,charsPer)
		 */

		template<typename T>
		void Encode( const Symbol s[], size_t len, size_t charsPerWord, vector<T> & code ) const {
			size_t size = symbols.size();
			size_t words = len - charsPerWord + 1;
			code.resize( words );

			for ( size_t offset = 0; offset < words; offset++ ) {
				code[offset] = Horner<T>( s + offset, size, charsPerWord );
			}
		}

		/**
		 *	Another alternative experimental coding which implements either a perfect kmer hash, 
		 *	if the result type is wide enough for the kmer, or a very bad hash otherwise.
		 *	Using the ambiguous amino acid alphabet with |\Sigma|=24, this will work well for 
		 *	1 <= len <= 13. With the ambiguous DNA alphabet |\Sigma|=8, so it will work well
		 *	for 1 <= len <= 21.
		 *
		 *	@typeparam	T Probably size_t, or something like that.
		 *	@param[in]	s A sequence of bytes to be converted into a sequence of packed n-grams.
		 *	@param[in]	len The number of symbols in s to process.
		 */

		template<typename T>
		T Encode( const Symbol s[], size_t len ) const {
			size_t size = symbols.size();
			return Horner<T>( s, size, len );
		}

		/**
		 *	Uses Horners method to pack a sequence of symbols into a word of some type.
		 *	@typeparam T The numeric type of desired result.
		 *	@param s A sequence of bytes to be turned into a number.
		 *	@param radix The radix of the numeric representation/encoding.
		 *	@param len The number of terms to pack into the number.
		 *	@returns s[0] \times radix^(l-1) + s[1] \times radix^(l-2) + ... + s[l-1]
		 */
		template<typename T>
		static T Horner( const Symbol s[], size_t radix, size_t len ) {
			T ret = 0;

			for ( size_t i = 0; i < len; i++ ) {
				ret = ret * radix + s[i].value;
			}

			return ret;
		}

		/// <summary> Encodes the designated string by computing a codeword for each tuple of 
		///		charsPerWord contiguous symbols. Multiple lists of codewords are calculated,
		///		one for each offset in 0, ..., (charsPerWord-1), and each of these lists becomes 
		///		a row of the resulting output dataset, code. 
		///		
		///		After processing:
		///		*	code[offset] contains the codewords of tuples starting at positions
		///			i*charsPerWord + offset for i = 1..(len / charsPerWord).
		/// </summary>
		/// <param name="s">The character sequence to encode.</param>
		/// <param name="charsPerWord">The number of encoded characters to pack into each word.
		///		This is independent of any particular kmer tiling scheme. It has to do with the
		///		maximum size of a pre-computed kmer similarity lookup table that can be created.
		///		<para>
		///			For my current experiments, I use charsPerWord \in {2,3} for amino acid, with the
		///			BLOSUM alphabet, this yields in-memory distance tables of approx 1/4GB.
		///		</para>
		///		<para>
		///			To encode DNA, I place one kmer per codeword, which introduces a constraint 
		///			on the maximum length, which becomes CHAR_BIT * sizeof(KmerWord) / 2.
		///		</para>
		/// </param>
		/// <param name="len">The number of characters to encode.</param>
		/// <param name="kmerLength">The number of characters in each kmer.</param>
		///	<param name="code">
		///		A reference to a list of arrays which will hold the codes 
		///		starting at each possible offset: 0 .. (charsPerWord - 1).
		///		<para>
		///			The staggered arrangement allows a sequence of successive kmers to be read from
		///			contiguous memory locations, so that code[0] is a list of partial kmers which start 
		///			at {0, charsPerWord, 2*charsPerWord, ...}; code[1] is the list of partial k-mers that
		///			start at {1, 1+charsPerWord, 1+2*charsPerWord, ...}.
		///		</para>
		///		<para>
		///			For DNA, where a KmerWord contains a whole k-mer, code will contain just a single
		///			row. That is, code[0] will contain the list of packed kmers.
		///		</para>
		///	</param>
		/// <returns></returns>

		virtual void Encode( 
			const Symbol s[], 
			size_t len, 
			size_t kmerLength, 
			size_t charsPerWord, 
			vector<vector<KmerWord>> & code 
		) const {
			// TODO: Get RID of this !!!
			// Assert::IsTrue( charsPerWord * this->BitsPerSymbol() <= KWORD_BITS, FileAndLine );

			if ( len + 1 < kmerLength ) {
				stringstream str;
				str << "Alphabet::Encode: string to encode must contain at least one k-mer: kmerLength = "
					<< kmerLength
					<< ", len = "
					<< len
					<< "\n";
				throw Exception( str.str(), FileAndLine );
			}

			if ( kmerLength > charsPerWord ) {
				if ( kmerLength % charsPerWord != 0 ) {
					stringstream str;
					str << "Alphabet::Encode: kmerLength must be divisible by charsPerWord: kmerLength = "
						<< kmerLength
						<< ", charsPerWord = "
						<< charsPerWord
						<< "\n"
						<< "The reason for this is that individual symbol bit patterns for the entire sequence are concatenated into one packed array for each offset starting at 0, 1, ..., (charsPerWord-1). If we try to use this with k-mer length not divisible by charsPerWord, we get garbage from adjacent kmers in the packed array."
						<< "\n";
					throw Exception( str.str(), FileAndLine );
				}

				code.resize( charsPerWord );

				for ( size_t i = 0; i < charsPerWord; i++ ) {
					size_t codeWordsNeeded = WordsRequiredToPack( len - i, charsPerWord );
					code[i].resize( codeWordsNeeded );
					Encode( s + i, len - i, charsPerWord, code[i].data() );
				}
			}
			else {
				// This is primarily for DNA;
				code.resize( 1 );
				code[0].clear();

				for ( size_t i = 0; i < len - kmerLength + 1; i++ ) {
					KmerWord codeWord;
					Encode( s + i, kmerLength, charsPerWord, &codeWord );
					code[0].push_back( codeWord );
				}
			}
		}

		/// <summary> Decodes a sequence of zero-origin numeric values into a string.
		/// </summary>
		/// <param name="code">The sequence of numeric code values.</param>
		/// <param name="k">The number of symbols required in the decoded string.</param>
		/// <param name="charsPerWord">The number of characters packed into each 
		///		word. See Encode for details.
		///	</param>
		/// <returns></returns>

		virtual void Decode( const KmerWord code[], size_t k, size_t charsPerWord, Symbol charBuffer[] ) const {
			// (cerr << "B").flush();
			size_t words = (k + charsPerWord - 1) / charsPerWord;
			size_t remaining = k % charsPerWord;

			if ( remaining == 0 ) {
				words++;
			}

			size_t size = symbols.size();
			Symbol* buffPtr = charBuffer;

			for ( size_t wordIndex = 1; wordIndex < words; wordIndex++, buffPtr += charsPerWord ) {
				KmerWord currentWord = code[wordIndex - 1];
				Symbol* last = buffPtr + charsPerWord - 1;

				for ( size_t i = 0; i < charsPerWord; i++, last-- ) {
					last->value = currentWord % size;
					currentWord /= size;
				}
			}

			if ( remaining > 0 ) {
				KmerWord currentWord = code[words - 1];
				Symbol* last = buffPtr + remaining - 1;

				for ( size_t i = 0; i < remaining; i++, last-- ) {
					last->value = currentWord % size;
					currentWord /= size;
				}
			}
		}

		/// <summary> Decodes a sequence of zero-origin numeric values into a string.
		/// </summary>
		/// <param name="code">The sequence of numeric code values.</param>
		/// <param name="len">The number of symbols required in the decoded string.</param>
		/// <param name="charsPerWord">The number of characters packed into each 
		///		word. See Pack for details.
		///	</param>
		/// <returns></returns>

		virtual void Decode( 
			const vector<vector<KmerWord>> & code, 
			size_t len, 
			size_t kmerLength, 
			size_t charsPerWord, 
			Symbol charBuffer[] 
		) const {
			// (cerr << "A").flush();
			if ( kmerLength > charsPerWord ) {
				Decode( code[0].data(), len, charsPerWord, charBuffer );
			}
			else {
				for ( size_t i = 0; i < len; i++ ) {
					Decode( &code[0][i], kmerLength, charsPerWord, charBuffer + i );
				}
			}
		}

		/**
		 *	<summary> Encodes the items into byte, on symbol per code byte.
		 *	The default encoding uses the ordinal position of each symbol in
		 *	Symbols() as its numeric code.
		 *	</summary>
		 *	<param name="s">A string to encode.</param>
		 *	<param name="code">a vector into which the encoded values will be saved.</param>
		 *	<returns></returns>
		 */

		virtual vector<Symbol> Encode( const string & s ) const {
			vector<Symbol> code( s.length() );

			for ( size_t i = 0; i < s.length(); i++ ) {
				code[i] = Encode( s[i] );
			}

			return code;
		}

		/**
		 *	<summary> Encodes the items into byte, on symbol per code byte.
		 *	The default encoding uses the ordinal position of each symbol in
		 *	Symbols() as its numeric code.
		 *	</summary>
		 *	<param name="s">A string to encode.</param>
		 *	<param name="code">a vector into which the encoded values will be saved.</param>
		 *	<returns></returns>
		 */

		virtual string Decode( const vector<Symbol> & code ) const {
			string s( code.size(), Decode(DefaultSymbol()) );

			for ( size_t i = 0; i < code.size(); i++ ) {
				s[i] = Decode( code[i] );
			}

			return s;
		}

		/**
		 *	<summary> Encodes the items into byte, on symbol per code byte.
		 *	The default encoding uses the ordinal position of each symbol in
		 *	Symbols() as its numeric code.
		 *	</summary>
		 *	<param name="s">A string to encode.</param>
		 *	<param name="code">a vector into which the encoded values will be saved.</param>
		 *	<returns></returns>
		 */

		virtual string Decode( SubVector<Symbol> & code ) const {
			string s( code.size(), Decode(DefaultSymbol()) );

			for ( size_t i = 0; i < code.size(); i++ ) {
				s[i] = Decode( code[i] );
			}

			return s;
		}

		/**
		 *	<summary> Get numeric code of a symbol.
		 *	</summary>
		 *	<param name="ch">The character to encode.</param>
		 *	<returns>inverse[(unsigned char) ch];
		 *	</returns>
		 */

		virtual Symbol Encode( char ch ) const {
			// (cerr << ch << " ... " << (int)(inverse[(unsigned char)ch]) << "\n").flush();
			return Symbol::From(inverse[(unsigned char) ch]);
		}

		/**
		 *	<summary>
		 *		Get symbol corresponding to a numeric code.
		 *	</summary>
		 *	<param name="code">
		 *		The numeric value to decode.
		 *	</param>
		 *	<returns>
		 *		symbols[code];
		 *	</returns>
		 */

		virtual char Decode( Symbol code ) const {
			return symbols[code.value];
		}

		/// <summary> Gets the number of symbols in this alphabet.
		/// </summary>

		int Size() const {
			return (int) symbols.size();
		}

		/// <summary> Gets the symbols defined in this alphabet.
		/// </summary>

		string Symbols() const {
			return symbols;
		}

		/// Get the name of the alphabet.
		string Name() const {
			return name;
		}

		/// <summary>Gets a symbol which can be used to pad strings.</summary>

		Symbol DefaultSymbol() const {
			return defaultSymbol;
		}

		const vector<byte> & Inverse() const {
			return inverse;
		}

		virtual ~Alphabet() {}

		virtual string ReverseComplement( const string & sequence ) const {
			return string( "" );
		}

		bool operator==( const Alphabet & other ) const {
			return symbols == other.symbols;
		}

		bool operator != ( const Alphabet & other ) const {
			return !operator==( other );
		}
	};

	class AA : public Alphabet {
		// const string AminoAcids = "arndcqeghilkmfpstwyvbzx*";
	public:
		AA() : Alphabet( "AA", "arndcqeghilkmfpstwyvbzx*" ) {
			DefaultSymbol( Encode('*') );
		}
	};

	class DefaultAlphabet : public Alphabet {
	public:
		DefaultAlphabet() : Alphabet( "DefaultAlphabet", Util::PrintableChars() ) {
			DefaultSymbol( Encode(' ') );
		}
	};

	class DNA : public Alphabet {

	public:
		DNA() : Alphabet( "DNA", "nacgt" ) {
			DefaultSymbol( Encode('n') );
		}

		string ReverseComplement( const string & sequence ) {
			int len = sequence.length();
			string complement( len + 1, 'n' );

			for ( int i = 0; i < len; i++ ) {
				char ch = sequence[i];
				char cc = ch == 'a' ? symbols.back() : ch == 'c' ? 'g' : ch == 'g' ? 'c' : ch == 'n' ? 'n' : 'a';
				complement[len - i - 1] = cc;
			}

			return complement;
		}

		void ReverseComplement( const byte source[], size_t len, byte dest[] ) {
			for ( size_t i = 0; i < len; i++ ) {
				byte ch = source[i];
				byte cc = ch == 1 ? 8 : ch == 2 ? 4 : ch == 4 ? 2 : ch == 0 ? 0 : 1;
				dest[len - i - 1] = cc;
			}
		}

	};

	class RNA : public DNA {
	public:
		RNA() {
			symbols[4] = 'u';
			inverse['U'] = 4;
			inverse['u'] = 4;
			inverse['T'] = 0;
			inverse['t'] = 0;
		}
	};

	class Alphabets {
	public:
		static Alphabet * AA() {
			static QutBio::AA value;
			return &value;
		}

		static Alphabet * DNA() {
			static QutBio::DNA value;
			return &value;
		}

		static Alphabet * RNA() {
			static QutBio::RNA value;
			return &value;
		}

		/**
		 * Gets the address of a static object representing the
		 */

		static Alphabet * DEFAULT() {
			static QutBio::DefaultAlphabet value;
			return &value;
		}

		/**
		 * Factory which gets a reference to an alphabet.
		 *
		 * @param alphaName_or_symbols A string which should be either one of the
		 * 		"standard" names {AA, DNA, RNA, DEFAULT }, or a literal list of
		 * 		symbols which will be used to define a custom alphabet.
		 *
		 * @returns The address of an Alphabet object. Custom objects are created
		 * 		afresh each time the function is called, using operator new, so
		 * 		you might want to consider deleting them. Standard alphabets are
		 * 		pointers to static objects, so you definitely do not want to delete
		 * 		a standard alphabet.
		 */

		static Alphabet * ByName( const string & alphaName_or_symbols ) {
			return alphaName_or_symbols == "AA" ? AA()
				: alphaName_or_symbols == "DNA" ? DNA()
				: alphaName_or_symbols == "RNA" ? RNA()
				: alphaName_or_symbols == "DefaultAlphabet" ? DEFAULT()
				: new Alphabet( "Custom", alphaName_or_symbols );
		}

		static map<string, pAlphabet> StandardAlphabets() {
			map<string, pAlphabet> result;
			result["AA"] = AA();
			result["DNA"] = DNA();
			result["RNA"] = RNA();
			result["DEFAULT"] = DEFAULT();
			return result;
		}
	};

	typedef Alphabet * pAlphabet;
}
