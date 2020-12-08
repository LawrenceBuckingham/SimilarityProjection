#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cmath>
#include <climits>
#include <fstream>
#include <cstdint>

#include "Alphabet.hpp"
#include "Assert.hpp"
#include "Delegates.hpp"
#include "DistanceType.hpp"
#include "Histogram.hpp"
#include "IntegerDistribution.hpp"
#include "Types.hpp"
#include "CsvIO.hpp"

namespace QutBio {

#define BAD_DIST numeric_limits<Distance>::min()
#define IS_BAD_DIST(x) (x == BAD_DIST) 
#define MAX_DIST numeric_limits<Distance>::max()

#define BAD_SIM numeric_limits<int>::min()
#define IS_BAD_SIM(x) (x == BAD_SIM) 

	typedef class SimilarityMatrix* pSimilarityMatrix;

	struct SimilarityMatrix :
		public virtual ICsvWriter,
		public virtual ICsvReader {
		QutBio::Alphabet* alphabet; //< Underlying symbol alphabet.
		int dict[128][128]; //< internal lookup table
		int maxValue = numeric_limits<int>::min(); //< Largest (most positive) similarity value observed in matrix.
		int minValue = numeric_limits<int>::max(); //< Smallest (least positive) similarity value observed in maatrix.
		string id;
		FlatMatrix<int> digrams;

		/**
		 *	Construct a similarity matrix to compare symbols from the given alphabet.
		 *	Alternatively, the alphabet can be inferred from a serialised similarity
		 *	lookup table.
		 *	@param	alphabet The (possibly null) address of an alphabet. If null, then
		 *					 the alphabet will be inferred during the Parse operation.
		 */
		SimilarityMatrix(
			QutBio::Alphabet* alphabet = 0
		) : alphabet( alphabet ) //
		{
			for ( int i = 0; i < 128; i++ ) {
				for ( int j = 0; j < 128; j++ ) {
					dict[i][j] = BAD_SIM;
				}
			}
		}

		void SetSimilarity( Symbol s, Symbol t, int value ) {
			assert_true( s.value < 128 && t.value < 128 );
			dict[s.value][t.value] = value;

			if ( maxValue < value ) {
				maxValue = value;
			}

			if ( minValue > value ) {
				minValue = value;
			}
		}

		Symbol LookupSymbol( char ch ) {
			if ( alphabet->Symbols().find( ch ) != string::npos ) {
				return alphabet->Encode( ch );
			}
			else if ( alphabet->Symbols().find( tolower( ch ) ) != string::npos ) {
				return alphabet->Encode( tolower( ch ) );
			}
			else {
				Symbol res{127};
				return res;
			}
		}

		SimilarityMatrix& Parse( istream& reader ) {
			maxValue = numeric_limits<int>::min();
			minValue = numeric_limits<int>::max();

			for ( uint i = 0; i < 128; i++ ) {
				for ( uint j = 0; j < 128; j++ ) {
					dict[i][j] = BAD_SIM;
				}
			}

			vector<Symbol> symbols;

			string currentLine;
			int rowCounter = 0;
			bool headerDone = false;

			while ( !( reader.bad() || reader.eof() ) ) {
				getline( reader, currentLine );
				currentLine = String::Trim( currentLine );

				if ( currentLine.size() == 0 ) break;

				if ( currentLine[0] == '#' ) continue;

				auto parts = String::Split( currentLine, ' ' );
				size_t len = parts.size();

				if ( !headerDone ) {
					headerDone = true;

					// process the heading row

					if ( alphabet ) {
						for ( uint i = 0; i < len; i++ ) {
							symbols.push_back( LookupSymbol( parts[i][0] ) );
						}
					}
					else {
						string letters;

						for ( uint i = 0; i < len; i++ ) {
							char ch = parts[i][0];
							letters.push_back( ch );
						}

						alphabet = new QutBio::Alphabet( "custom", letters );
						symbols = alphabet->Encode( letters );
					}
				}
				else {
					int firstIndex = 0;
					Symbol rowSymbol = LookupSymbol( parts[0][0] );

					if ( rowSymbol.value != 127 ) {
						// Each row has a label: use it to specify first element of pair.
						firstIndex++;
					}
					else {
						// Rows contain no labels; assume same order as symbols appear in heading 
						// row.
						rowSymbol = symbols[rowCounter];
						rowCounter++;
					}

					// process numeric entries in row.
					for ( uint j = firstIndex; j < len; j++ ) {
						int d = (int) Double::Parse( parts[j] );
						SetSimilarity( rowSymbol, symbols[j - firstIndex], d );
					}
				}
			}

			for ( int i = 0; i < 128; i++ ) {
				for ( int j = 0; j < 128; j++ ) {
					if ( dict[i][j] == BAD_SIM ) {
						dict[i][j] = minValue;
					}
				}
			}

			ComputeDigramDifference( symbols );

			return *this;
		}

		void ComputeDigramDifference( const vector<Symbol>& symbols ) {
			uint8_t n = symbols.size();
			uint16_t digramVocabSize = n * n;
			digrams.resize( digramVocabSize, digramVocabSize );
			Symbol x[2], y[2];

			for ( x[0].value = 0; x[0].value < n; x[0].value++ ) {
				for ( x[1].value = 0; x[1].value < n; x[1].value++ ) {
					uint xOff = Alphabet::Horner<uint>( x, n, 2 );

					for ( y[0].value = 0; y[0].value < n; y[0].value++ ) {
						for ( y[1].value = 0; y[1].value < n; y[1].value++ ) {
							uint yOff = Alphabet::Horner<uint>( y, n, 2 );

							digrams( xOff, yOff ) = Difference( x, y, 2 );
						}
					}

				}
			}
		}

		/**
		 *	Serialise the similarity matrix to a CSV stream, fulfilling the ICsvWriter contract.
		 *	Equivalent to w.Write(ToString()).
		 *	@param w A CsvWriter which will serialise the content of this object.
		 */
		void Write( CsvWriter& w ) const override {
			w.Write( ToString() );
		}

		/**
		 *	Deserialise the similarity matrix from a CSV stream, fulfilling the ICsvReader contract.
		 *	@param r reference to a CsvReader from which data will be obtained.
		 */
		void Read( CsvReader& r ) override {
			string s;
			r >> s;
			istringstream str( s );
			Parse( str );
		}

		/**
		 *	Returns a string which can be parsed later to recreate this similarity matrix.
		 */
		string ToString() const {
			string s = alphabet->Symbols();
			ostringstream f;

			for ( auto p = s.begin(); p != s.end(); p++ ) {
				f << " " << *p;
			}

			for ( auto q = s.begin(); q < s.end(); q++ ) {
				f << "\n";

				for ( auto p = s.begin(); p != s.end(); p++ ) {
					f << " " << ( (int) Similarity( alphabet->Encode( *q ), alphabet->Encode( *p ) ) );
				}
			}

			return f.str();
		}

		SimilarityMatrix( const SimilarityMatrix& other ) {
			memcpy( dict, other.dict, sizeof( dict ) );
			maxValue = other.maxValue;
			minValue = other.minValue;
			id = other.id;
			digrams = other.digrams;
		}

		bool operator!=( const SimilarityMatrix& other ) const {
			return !operator==( other );
		}

		bool operator==( const SimilarityMatrix& other ) const {
			if ( *alphabet != *other.alphabet ) return false;

			if ( memcmp( dict, other.dict, sizeof( dict ) ) != 0 ) return false;

			return true;
		}

		SimilarityMatrix& operator=( const SimilarityMatrix& other ) {
			memcpy( dict, other.dict, sizeof( dict ) );
			maxValue = other.maxValue;
			minValue = other.minValue;
			id = other.id;
			digrams = other.digrams;
			return *this;
		}

		int8_t Similarity( Symbol s, Symbol t ) const {
			return dict[s.value][t.value];
		}

		int8_t MaxValue( void ) const {
			return maxValue;
		}

		int8_t MinValue( void ) const {
			return minValue;
		}

		void Foreach( Action3<int, int, Distance> action ) const {
			for ( int i = 0; i < 128; i++ ) {
				for ( int j = 0; j < 128; j++ ) {
					if ( !IS_BAD_SIM( dict[i][j] ) ) {
						action( i, j, dict[i][j] );
					}
				}
			}
		}

		Distance Similarity( const Symbol x[], const Symbol y[], uint length ) const {
			Distance score = 0;

			for ( uint i = 0; i < length; i++ ) {
				score += dict[x[i].value][y[i].value];
			}

			return score;
		}

		Distance Similarity( const char x[], const char y[], uint length ) const {
			Distance score = 0;

			for ( uint i = 0; i < length; i++ ) {
				score += dict[alphabet->Encode( x[i] ).value][alphabet->Encode( y[i] ).value];
			}

			return score;
		}

		Distance Similarity( const Symbol x[], uint length ) const {
			Distance score = 0;

			for ( uint i = 0; i < length; i++ ) {
				auto t = x[i].value;
				score += dict[t][t];
			}

			return score;
		}

		Distance HalperinDistance( const Symbol x[], const Symbol y[], int length ) const {
			Distance d = Similarity( x, length );
			d += Similarity( y, length );
			d -= 2 * Similarity( x, y, length );
			return d;
		}

		Distance Difference( const char x[], const char y[], uint length ) const {
			Distance d = length * maxValue - Similarity( x, y, length );
			return d;
		}

		Distance Difference( const Symbol x[], const Symbol y[], uint length ) const {
			Distance d = length * maxValue - Similarity( x, y, length );
			return d;
		}

		Distance Difference( Symbol x, Symbol y ) const {
			Distance d = maxValue - Similarity( x, y );
			return d;
		}

		Distance DigramDifference( const Digram* x, const Digram* y, uint len ) const {
			Distance ret = digrams.at( x[0], y[0] );

			for ( uint i = 1; i < len; i++ ) {
				uint offset = 2 * i;
				ret += digrams.at( x[offset], y[offset] );
			}

			return ret;
		}

		Distance DigramDifference( const Digram x, const Digram y ) const {
			Distance ret = digrams.at( x, y );
			return ret;
		}

		bool IsWithin( const byte* x, const byte* y, int length, Distance threshold, Distance& result ) {
			Distance dist = 0;

			for ( int i = 0; i < length; i++ ) {
				dist += maxValue - dict[x[i]][y[i]];

				if ( dist > threshold ) {
					return false;
				}
			}

			result = dist;
			return true;
		}

		static SimilarityMatrix* Blosum100() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum100.iij							\n"
				"#  * column uses minimum score 										\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 100											\n"
				"#  Entropy =   1.4516, Expected =  -1.0948								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -2 -3 -1 -1 -2 -1 -3 -3 -3 -2 -2 -4 -1  1 -1 -4 -4 -1 -3 -2 -1 -7\n"
				"-2  7 -1 -3 -5  0 -2 -4 -1 -4 -4  2 -2 -4 -3 -2 -2 -4 -3 -4 -2 -1 -2 -7\n"
				"-2 -1  7  1 -4 -1 -1 -2  0 -5 -5 -1 -4 -5 -4  0 -1 -6 -3 -4  4 -1 -2 -7\n"
				"-3 -3  1  7 -5 -2  1 -3 -2 -6 -6 -2 -5 -5 -3 -1 -2 -7 -5 -5  4  0 -3 -7\n"
				"-1 -5 -4 -5  9 -5 -6 -5 -5 -2 -3 -5 -3 -3 -5 -2 -2 -5 -4 -2 -5 -6 -3 -7\n"
				"-1  0 -1 -2 -5  7  1 -3  0 -4 -3  1 -1 -4 -2 -1 -2 -3 -3 -3 -1  3 -2 -7\n"
				"-2 -2 -1  1 -6  1  6 -4 -1 -5 -5  0 -4 -5 -3 -1 -2 -5 -4 -3  0  5 -2 -7\n"
				"-1 -4 -2 -3 -5 -3 -4  6 -4 -6 -5 -3 -5 -5 -4 -1 -3 -5 -6 -5 -2 -4 -3 -7\n"
				"-3 -1  0 -2 -5  0 -1 -4  9 -5 -4 -2 -3 -2 -3 -2 -3 -3  1 -5 -1 -1 -2 -7\n"
				"-3 -4 -5 -6 -2 -4 -5 -6 -5  5  1 -4  1 -1 -4 -4 -2 -4 -3  2 -5 -4 -2 -7\n"
				"-3 -4 -5 -6 -3 -3 -5 -5 -4  1  5 -4  2  0 -4 -4 -3 -4 -3  0 -5 -4 -2 -7\n"
				"-2  2 -1 -2 -5  1  0 -3 -2 -4 -4  6 -2 -4 -2 -1 -2 -5 -4 -4 -1  0 -2 -7\n"
				"-2 -2 -4 -5 -3 -1 -4 -5 -3  1  2 -2  8 -1 -4 -3 -2 -3 -3  0 -4 -3 -2 -7\n"
				"-4 -4 -5 -5 -3 -4 -5 -5 -2 -1  0 -4 -1  7 -5 -3 -3  0  3 -2 -5 -5 -3 -7\n"
				"-1 -3 -4 -3 -5 -2 -3 -4 -3 -4 -4 -2 -4 -5  8 -2 -3 -6 -5 -4 -3 -3 -3 -7\n"
				" 1 -2  0 -1 -2 -1 -1 -1 -2 -4 -4 -1 -3 -3 -2  6  1 -4 -3 -3 -1 -1 -1 -7\n"
				"-1 -2 -1 -2 -2 -2 -2 -3 -3 -2 -3 -2 -2 -3 -3  1  6 -5 -3 -1 -2 -2 -1 -7\n"
				"-4 -4 -6 -7 -5 -3 -5 -5 -3 -4 -4 -5 -3  0 -6 -4 -5 11  1 -4 -6 -4 -4 -7\n"
				"-4 -3 -3 -5 -4 -3 -4 -6  1 -3 -3 -4 -3  3 -5 -3 -3  1  8 -3 -4 -4 -3 -7\n"
				"-1 -4 -4 -5 -2 -3 -3 -5 -5  2  0 -4  0 -2 -4 -3 -1 -4 -3  5 -5 -3 -2 -7\n"
				"-3 -2  4  4 -5 -1  0 -2 -1 -5 -5 -1 -4 -5 -3 -1 -2 -6 -4 -5  4  0 -2 -7\n"
				"-2 -1 -1  0 -6  3  5 -4 -1 -4 -4  0 -3 -5 -3 -1 -2 -4 -4 -3  0  4 -2 -7\n"
				"-1 -2 -2 -3 -3 -2 -2 -3 -2 -2 -2 -2 -2 -3 -3 -1 -1 -4 -3 -2 -2 -2 -2 -7\n"
				"-7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::AA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		static SimilarityMatrix* Blosum80() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum80.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 80											\n"
				"#  Entropy =   0.9868, Expected =  -0.7442								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -2 -2 -1 -1 -1  0 -2 -2 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1 -6\n"
				"-2  6 -1 -2 -4  1 -1 -3  0 -3 -3  2 -2 -4 -2 -1 -1 -4 -3 -3 -2  0 -1 -6\n"
				"-2 -1  6  1 -3  0 -1 -1  0 -4 -4  0 -3 -4 -3  0  0 -4 -3 -4  4  0 -1 -6\n"
				"-2 -2  1  6 -4 -1  1 -2 -2 -4 -5 -1 -4 -4 -2 -1 -1 -6 -4 -4  4  1 -2 -6\n"
				"-1 -4 -3 -4  9 -4 -5 -4 -4 -2 -2 -4 -2 -3 -4 -2 -1 -3 -3 -1 -4 -4 -3 -6\n"
				"-1  1  0 -1 -4  6  2 -2  1 -3 -3  1  0 -4 -2  0 -1 -3 -2 -3  0  3 -1 -6\n"
				"-1 -1 -1  1 -5  2  6 -3  0 -4 -4  1 -2 -4 -2  0 -1 -4 -3 -3  1  4 -1 -6\n"
				" 0 -3 -1 -2 -4 -2 -3  6 -3 -5 -4 -2 -4 -4 -3 -1 -2 -4 -4 -4 -1 -3 -2 -6\n"
				"-2  0  0 -2 -4  1  0 -3  8 -4 -3 -1 -2 -2 -3 -1 -2 -3  2 -4 -1  0 -2 -6\n"
				"-2 -3 -4 -4 -2 -3 -4 -5 -4  5  1 -3  1 -1 -4 -3 -1 -3 -2  3 -4 -4 -2 -6\n"
				"-2 -3 -4 -5 -2 -3 -4 -4 -3  1  4 -3  2  0 -3 -3 -2 -2 -2  1 -4 -3 -2 -6\n"
				"-1  2  0 -1 -4  1  1 -2 -1 -3 -3  5 -2 -4 -1 -1 -1 -4 -3 -3 -1  1 -1 -6\n"
				"-1 -2 -3 -4 -2  0 -2 -4 -2  1  2 -2  6  0 -3 -2 -1 -2 -2  1 -3 -2 -1 -6\n"
				"-3 -4 -4 -4 -3 -4 -4 -4 -2 -1  0 -4  0  6 -4 -3 -2  0  3 -1 -4 -4 -2 -6\n"
				"-1 -2 -3 -2 -4 -2 -2 -3 -3 -4 -3 -1 -3 -4  8 -1 -2 -5 -4 -3 -2 -2 -2 -6\n"
				" 1 -1  0 -1 -2  0  0 -1 -1 -3 -3 -1 -2 -3 -1  5  1 -4 -2 -2  0  0 -1 -6\n"
				" 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -2 -1 -1 -2 -2  1  5 -4 -2  0 -1 -1 -1 -6\n"
				"-3 -4 -4 -6 -3 -3 -4 -4 -3 -3 -2 -4 -2  0 -5 -4 -4 11  2 -3 -5 -4 -3 -6\n"
				"-2 -3 -3 -4 -3 -2 -3 -4  2 -2 -2 -3 -2  3 -4 -2 -2  2  7 -2 -3 -3 -2 -6\n"
				" 0 -3 -4 -4 -1 -3 -3 -4 -4  3  1 -3  1 -1 -3 -2  0 -3 -2  4 -4 -3 -1 -6\n"
				"-2 -2  4  4 -4  0  1 -1 -1 -4 -4 -1 -3 -4 -2  0 -1 -5 -3 -4  4  0 -2 -6\n"
				"-1  0  0  1 -4  3  4 -3  0 -4 -3  1 -2 -4 -2  0 -1 -4 -3 -3  0  4 -1 -6\n"
				"-1 -1 -1 -2 -3 -1 -1 -2 -2 -2 -2 -1 -1 -2 -2 -1 -1 -3 -2 -1 -2 -1 -1 -6\n"
				"-6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::AA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		static SimilarityMatrix* Blosum62() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum62.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 62											\n"
				"#  Entropy =   0.6979, Expected =  -0.5209								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n"
				"-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n"
				"-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n"
				"-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n"
				" 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n"
				"-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n"
				"-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n"
				" 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n"
				"-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n"
				"-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n"
				"-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n"
				"-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n"
				"-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n"
				"-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n"
				"-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n"
				" 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n"
				" 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n"
				"-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n"
				"-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n"
				" 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n"
				"-2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  0 -1 -4\n"
				"-1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  0  4 -1 -4\n"
				" 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n"
				"-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::AA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		static SimilarityMatrix* Blosum50() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum50.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 50											\n"
				"#  Entropy =   0.4808, Expected =  -0.3573								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -1 -2 -1 -1 -1  0 -2 -1 -2 -1 -1 -3 -1  1  0 -3 -2  0 -2 -1 -1 -5\n"
				"-2  7 -1 -2 -4  1  0 -3  0 -4 -3  3 -2 -3 -3 -1 -1 -3 -1 -3 -1  0 -1 -5\n"
				"-1 -1  7  2 -2  0  0  0  1 -3 -4  0 -2 -4 -2  1  0 -4 -2 -3  4  0 -1 -5\n"
				"-2 -2  2  8 -4  0  2 -1 -1 -4 -4 -1 -4 -5 -1  0 -1 -5 -3 -4  5  1 -1 -5\n"
				"-1 -4 -2 -4 13 -3 -3 -3 -3 -2 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -3 -3 -2 -5\n"
				"-1  1  0  0 -3  7  2 -2  1 -3 -2  2  0 -4 -1  0 -1 -1 -1 -3  0  4 -1 -5\n"
				"-1  0  0  2 -3  2  6 -3  0 -4 -3  1 -2 -3 -1 -1 -1 -3 -2 -3  1  5 -1 -5\n"
				" 0 -3  0 -1 -3 -2 -3  8 -2 -4 -4 -2 -3 -4 -2  0 -2 -3 -3 -4 -1 -2 -2 -5\n"
				"-2  0  1 -1 -3  1  0 -2 10 -4 -3  0 -1 -1 -2 -1 -2 -3  2 -4  0  0 -1 -5\n"
				"-1 -4 -3 -4 -2 -3 -4 -4 -4  5  2 -3  2  0 -3 -3 -1 -3 -1  4 -4 -3 -1 -5\n"
				"-2 -3 -4 -4 -2 -2 -3 -4 -3  2  5 -3  3  1 -4 -3 -1 -2 -1  1 -4 -3 -1 -5\n"
				"-1  3  0 -1 -3  2  1 -2  0 -3 -3  6 -2 -4 -1  0 -1 -3 -2 -3  0  1 -1 -5\n"
				"-1 -2 -2 -4 -2  0 -2 -3 -1  2  3 -2  7  0 -3 -2 -1 -1  0  1 -3 -1 -1 -5\n"
				"-3 -3 -4 -5 -2 -4 -3 -4 -1  0  1 -4  0  8 -4 -3 -2  1  4 -1 -4 -4 -2 -5\n"
				"-1 -3 -2 -1 -4 -1 -1 -2 -2 -3 -4 -1 -3 -4 10 -1 -1 -4 -3 -3 -2 -1 -2 -5\n"
				" 1 -1  1  0 -1  0 -1  0 -1 -3 -3  0 -2 -3 -1  5  2 -4 -2 -2  0  0 -1 -5\n"
				" 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  2  5 -3 -2  0  0 -1  0 -5\n"
				"-3 -3 -4 -5 -5 -1 -3 -3 -3 -3 -2 -3 -1  1 -4 -4 -3 15  2 -3 -5 -2 -3 -5\n"
				"-2 -1 -2 -3 -3 -1 -2 -3  2 -1 -1 -2  0  4 -3 -2 -2  2  8 -1 -3 -2 -1 -5\n"
				" 0 -3 -3 -4 -1 -3 -3 -4 -4  4  1 -3  1 -1 -3 -2  0 -3 -1  5 -4 -3 -1 -5\n"
				"-2 -1  4  5 -3  0  1 -1  0 -4 -4  0 -3 -4 -2  0  0 -5 -3 -4  5  0 -1 -5\n"
				"-1  0  0  1 -3  4  5 -2  0 -3 -3  1 -1 -4 -1  0 -1 -2 -2 -3  0  5 -1 -5\n"
				"-1 -1 -1 -1 -2 -1 -1 -2 -1 -1 -1 -1 -1 -2 -2 -1  0 -3 -1 -1 -1 -1 -1 -5\n"
				"-5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::AA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		static SimilarityMatrix* Blosum45() {
#pragma region Data
			static const string data = "#  Matrix made by matblas from blosum45.iij\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/3 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 45											\n"
				"#  Entropy =   0.3795, Expected =  -0.2789								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -1 -2 -1 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -2 -2  0 -1 -1  0 -5\n"
				"-2  7  0 -1 -3  1  0 -2  0 -3 -2  3 -1 -2 -2 -1 -1 -2 -1 -2 -1  0 -1 -5\n"
				"-1  0  6  2 -2  0  0  0  1 -2 -3  0 -2 -2 -2  1  0 -4 -2 -3  4  0 -1 -5\n"
				"-2 -1  2  7 -3  0  2 -1  0 -4 -3  0 -3 -4 -1  0 -1 -4 -2 -3  5  1 -1 -5\n"
				"-1 -3 -2 -3 12 -3 -3 -3 -3 -3 -2 -3 -2 -2 -4 -1 -1 -5 -3 -1 -2 -3 -2 -5\n"
				"-1  1  0  0 -3  6  2 -2  1 -2 -2  1  0 -4 -1  0 -1 -2 -1 -3  0  4 -1 -5\n"
				"-1  0  0  2 -3  2  6 -2  0 -3 -2  1 -2 -3  0  0 -1 -3 -2 -3  1  4 -1 -5\n"
				" 0 -2  0 -1 -3 -2 -2  7 -2 -4 -3 -2 -2 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -5\n"
				"-2  0  1  0 -3  1  0 -2 10 -3 -2 -1  0 -2 -2 -1 -2 -3  2 -3  0  0 -1 -5\n"
				"-1 -3 -2 -4 -3 -2 -3 -4 -3  5  2 -3  2  0 -2 -2 -1 -2  0  3 -3 -3 -1 -5\n"
				"-1 -2 -3 -3 -2 -2 -2 -3 -2  2  5 -3  2  1 -3 -3 -1 -2  0  1 -3 -2 -1 -5\n"
				"-1  3  0  0 -3  1  1 -2 -1 -3 -3  5 -1 -3 -1 -1 -1 -2 -1 -2  0  1 -1 -5\n"
				"-1 -1 -2 -3 -2  0 -2 -2  0  2  2 -1  6  0 -2 -2 -1 -2  0  1 -2 -1 -1 -5\n"
				"-2 -2 -2 -4 -2 -4 -3 -3 -2  0  1 -3  0  8 -3 -2 -1  1  3  0 -3 -3 -1 -5\n"
				"-1 -2 -2 -1 -4 -1  0 -2 -2 -2 -3 -1 -2 -3  9 -1 -1 -3 -3 -3 -2 -1 -1 -5\n"
				" 1 -1  1  0 -1  0  0  0 -1 -2 -3 -1 -2 -2 -1  4  2 -4 -2 -1  0  0  0 -5\n"
				" 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -1 -1  2  5 -3 -1  0  0 -1  0 -5\n"
				"-2 -2 -4 -4 -5 -2 -3 -2 -3 -2 -2 -2 -2  1 -3 -4 -3 15  3 -3 -4 -2 -2 -5\n"
				"-2 -1 -2 -2 -3 -1 -2 -3  2  0  0 -1  0  3 -3 -2 -1  3  8 -1 -2 -2 -1 -5\n"
				" 0 -2 -3 -3 -1 -3 -3 -3 -3  3  1 -2  1  0 -3 -1  0 -3 -1  5 -3 -3 -1 -5\n"
				"-1 -1  4  5 -2  0  1 -1  0 -3 -3  0 -2 -3 -2  0  0 -4 -2 -3  4  0 -1 -5\n"
				"-1  0  0  1 -3  4  4 -2  0 -3 -2  1 -1 -3 -1  0 -1 -2 -2 -3  0  4 -1 -5\n"
				" 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1  0  0 -2 -1 -1 -1 -1 -1 -5\n"
				"-5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::AA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		static SimilarityMatrix* Blosum40() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum40.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/4 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 40											\n"
				"#  Entropy =   0.2851, Expected =  -0.2090								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -2 -1 -1 -2  0 -1  1 -2 -1 -2 -1 -1 -3 -2  1  0 -3 -2  0 -1 -1  0 -6\n"
				"-2  9  0 -1 -3  2 -1 -3  0 -3 -2  3 -1 -2 -3 -1 -2 -2 -1 -2 -1  0 -1 -6\n"
				"-1  0  8  2 -2  1 -1  0  1 -2 -3  0 -2 -3 -2  1  0 -4 -2 -3  4  0 -1 -6\n"
				"-1 -1  2  9 -2 -1  2 -2  0 -4 -3  0 -3 -4 -2  0 -1 -5 -3 -3  6  1 -1 -6\n"
				"-2 -3 -2 -2 16 -4 -2 -3 -4 -4 -2 -3 -3 -2 -5 -1 -1 -6 -4 -2 -2 -3 -2 -6\n"
				" 0  2  1 -1 -4  8  2 -2  0 -3 -2  1 -1 -4 -2  1 -1 -1 -1 -3  0  4 -1 -6\n"
				"-1 -1 -1  2 -2  2  7 -3  0 -4 -2  1 -2 -3  0  0 -1 -2 -2 -3  1  5 -1 -6\n"
				" 1 -3  0 -2 -3 -2 -3  8 -2 -4 -4 -2 -2 -3 -1  0 -2 -2 -3 -4 -1 -2 -1 -6\n"
				"-2  0  1  0 -4  0  0 -2 13 -3 -2 -1  1 -2 -2 -1 -2 -5  2 -4  0  0 -1 -6\n"
				"-1 -3 -2 -4 -4 -3 -4 -4 -3  6  2 -3  1  1 -2 -2 -1 -3  0  4 -3 -4 -1 -6\n"
				"-2 -2 -3 -3 -2 -2 -2 -4 -2  2  6 -2  3  2 -4 -3 -1 -1  0  2 -3 -2 -1 -6\n"
				"-1  3  0  0 -3  1  1 -2 -1 -3 -2  6 -1 -3 -1  0  0 -2 -1 -2  0  1 -1 -6\n"
				"-1 -1 -2 -3 -3 -1 -2 -2  1  1  3 -1  7  0 -2 -2 -1 -2  1  1 -3 -2  0 -6\n"
				"-3 -2 -3 -4 -2 -4 -3 -3 -2  1  2 -3  0  9 -4 -2 -1  1  4  0 -3 -4 -1 -6\n"
				"-2 -3 -2 -2 -5 -2  0 -1 -2 -2 -4 -1 -2 -4 11 -1  0 -4 -3 -3 -2 -1 -2 -6\n"
				" 1 -1  1  0 -1  1  0  0 -1 -2 -3  0 -2 -2 -1  5  2 -5 -2 -1  0  0  0 -6\n"
				" 0 -2  0 -1 -1 -1 -1 -2 -2 -1 -1  0 -1 -1  0  2  6 -4 -1  1  0 -1  0 -6\n"
				"-3 -2 -4 -5 -6 -1 -2 -2 -5 -3 -1 -2 -2  1 -4 -5 -4 19  3 -3 -4 -2 -2 -6\n"
				"-2 -1 -2 -3 -4 -1 -2 -3  2  0  0 -1  1  4 -3 -2 -1  3  9 -1 -3 -2 -1 -6\n"
				" 0 -2 -3 -3 -2 -3 -3 -4 -4  4  2 -2  1  0 -3 -1  1 -3 -1  5 -3 -3 -1 -6\n"
				"-1 -1  4  6 -2  0  1 -1  0 -3 -3  0 -3 -3 -2  0  0 -4 -3 -3  5  0 -1 -6\n"
				"-1  0  0  1 -3  4  5 -2  0 -4 -2  1 -2 -4 -1  0 -1 -2 -2 -3  0  5 -1 -6\n"
				" 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1  0 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -6\n"
				"-6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6 -6  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::AA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		static SimilarityMatrix* Blosum35() {
#pragma region Data
			static const string data =
				"#  Matrix made by matblas from blosum35.iij							\n"
				"#  * column uses minimum score											\n"
				"#  BLOSUM Clustered Scoring Matrix in 1/4 Bit Units					\n"
				"#  Blocks Database = /data/blocks_5.0/blocks.dat						\n"
				"#  Cluster Percentage: >= 35											\n"
				"#  Entropy =   0.2111, Expected =  -0.1550								\n"
				" A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\n"
				" 5 -1 -1 -1 -2  0 -1  0 -2 -1 -2  0  0 -2 -2  1  0 -2 -1  0 -1 -1  0 -5\n"
				"-1  8 -1 -1 -3  2 -1 -2 -1 -3 -2  2  0 -1 -2 -1 -2  0  0 -1 -1  0 -1 -5\n"
				"-1 -1  7  1 -1  1 -1  1  1 -1 -2  0 -1 -1 -2  0  0 -2 -2 -2  4  0  0 -5\n"
				"-1 -1  1  8 -3 -1  2 -2  0 -3 -2 -1 -3 -3 -1 -1 -1 -3 -2 -2  5  1 -1 -5\n"
				"-2 -3 -1 -3 15 -3 -1 -3 -4 -4 -2 -2 -4 -4 -4 -3 -1 -5 -5 -2 -2 -2 -2 -5\n"
				" 0  2  1 -1 -3  7  2 -2 -1 -2 -2  0 -1 -4  0  0  0 -1  0 -3  0  4 -1 -5\n"
				"-1 -1 -1  2 -1  2  6 -2 -1 -3 -1  1 -2 -3  0  0 -1 -1 -1 -2  0  5 -1 -5\n"
				" 0 -2  1 -2 -3 -2 -2  7 -2 -3 -3 -1 -1 -3 -2  1 -2 -1 -2 -3  0 -2 -1 -5\n"
				"-2 -1  1  0 -4 -1 -1 -2 12 -3 -2 -2  1 -3 -1 -1 -2 -4  0 -4  0 -1 -1 -5\n"
				"-1 -3 -1 -3 -4 -2 -3 -3 -3  5  2 -2  1  1 -1 -2 -1 -1  0  4 -2 -3  0 -5\n"
				"-2 -2 -2 -2 -2 -2 -1 -3 -2  2  5 -2  3  2 -3 -2  0  0  0  2 -2 -2  0 -5\n"
				" 0  2  0 -1 -2  0  1 -1 -2 -2 -2  5  0 -1  0  0  0  0 -1 -2  0  1  0 -5\n"
				" 0  0 -1 -3 -4 -1 -2 -1  1  1  3  0  6  0 -3 -1  0  1  0  1 -2 -2  0 -5\n"
				"-2 -1 -1 -3 -4 -4 -3 -3 -3  1  2 -1  0  8 -4 -1 -1  1  3  1 -2 -3 -1 -5\n"
				"-2 -2 -2 -1 -4  0  0 -2 -1 -1 -3  0 -3 -4 10 -2  0 -4 -3 -3 -1  0 -1 -5\n"
				" 1 -1  0 -1 -3  0  0  1 -1 -2 -2  0 -1 -1 -2  4  2 -2 -1 -1  0  0  0 -5\n"
				" 0 -2  0 -1 -1  0 -1 -2 -2 -1  0  0  0 -1  0  2  5 -2 -2  1 -1 -1  0 -5\n"
				"-2  0 -2 -3 -5 -1 -1 -1 -4 -1  0  0  1  1 -4 -2 -2 16  3 -2 -3 -1 -1 -5\n"
				"-1  0 -2 -2 -5  0 -1 -2  0  0  0 -1  0  3 -3 -1 -2  3  8  0 -2 -1 -1 -5\n"
				" 0 -1 -2 -2 -2 -3 -2 -3 -4  4  2 -2  1  1 -3 -1  1 -2  0  5 -2 -2  0 -5\n"
				"-1 -1  4  5 -2  0  0  0  0 -2 -2  0 -2 -2 -1  0 -1 -3 -2 -2  5  0 -1 -5\n"
				"-1  0  0  1 -2  4  5 -2 -1 -3 -2  1 -2 -3  0  0 -1 -1 -1 -2  0  4  0 -5\n"
				" 0 -1  0 -1 -2 -1 -1 -1 -1  0  0  0  0 -1 -1  0  0 -1 -1  0 -1  0 -1 -5\n"
				"-5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5 -5  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::AA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		static SimilarityMatrix* GetBlosum( int matrixId ) {
			switch ( matrixId ) {
			case 100: return Blosum100();
			case 80: return Blosum80();
			case 62: return Blosum62();
			case 50: return Blosum50();
			case 45: return Blosum45();
			case 40: return Blosum40();
			case 35: return Blosum35();
			default: return 0;
			}
		}

		/**
		 *	<summary>
		 *		Computes the value of the Karlin-Altshul parameter lambda for stipulated
		 *		background distributions.
		 *	</summary>
		 *	<param name=p1>
		 *		A normalised histogram  which contains the distribution of query residues.
		 *	</summary>
		 *	<param name=p2>
		 *		A normalised histogram  which contains the distribution of subject residues.
		 *	</summary>
		 */
		double ComputeExtremeDistLambda(
			const Histogram<char>& p1,
			const Histogram<char>& p2
		) {
			auto f = [=]( double lambda ) {
				double sum = 0;

				for ( auto x : p1.data ) {
					for ( auto y : p2.data ) {
						double score = dict[(uint8_t) x.first][(uint8_t) y.first];
						sum += x.second * y.second * exp( lambda * score );
					}
				}

				return sum;
			};

			auto fPrime = [=]( double lambda ) {
				double sum = 0;

				for ( auto x : p1.data ) {
					for ( auto y : p2.data ) {
						double score = dict[(uint8_t) x.first][(uint8_t) y.first];
						sum += x.second * y.second * score * exp( lambda * score );
					}
				}

				return sum;
			};

			double lambda = 0;

			while ( true ) {
				double lambda2 = lambda - f( lambda ) / fPrime( lambda );
				const double eps = 1e-10;

				if ( fabs( lambda2 - lambda ) < eps ) break;

				lambda = lambda2;
			}

			return lambda;
		}

		static SimilarityMatrix* GetDnaMatrix() {
#pragma region Data
			static const string data =
				"#  Matrix invented by Lawrence Buckingham. \n"
				" A  C  G  T  N\n"
				"16  0  0  0  4\n"
				" 0 16  0  0  4\n"
				" 0  0 16  0  4\n"
				" 0  0  0 16  4\n"
				" 4  4  4  4  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::DNA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		static SimilarityMatrix* GetRnaMatrix() {
#pragma region Data
			static const string data =
				"#  Matrix invented by Lawrence Buckingham. \n"
				" A  C  G  U  N\n"
				"16  0  0  0  4\n"
				" 0 16  0  0  4\n"
				" 0  0 16  0  4\n"
				" 0  0  0 16  4\n"
				" 4  4  4  4  1";
#pragma endregion
			static SimilarityMatrix matrix( Alphabets::RNA() );

			if ( IS_BAD_SIM( matrix.dict['a']['a'] ) ) {
				istringstream reader( data );
				matrix.Parse( reader );
			}

			return &matrix;
		}

		/**
		 *	<summary>
		 *		Gets either an in-built or custom similarity matrix defined by the supplied arguments.
		 *
		 *		Beware: if Custom is used, a new instance is created via the new operator. You will have to free this yourself when finished with it.
		 *	</summary>
		 */
		static SimilarityMatrix* GetMatrix(
			Alphabet* alphabet,
			DistanceType* dist,
			int id,
			const string& customFileName
		) {
			if ( dist == DistanceType::HalperinEtAl() || dist == DistanceType::BlosumDistance() ) {
				return GetBlosum( id );
			}
			else if ( dist == DistanceType::Custom() ) {
				SimilarityMatrix* matrix = new SimilarityMatrix( alphabet );

				if ( matrix == 0 ) {
					throw Exception( "Out of memory!", FileAndLine );
				}

				ifstream f( customFileName );
				matrix->Parse( f );
				return matrix;
			}
			else {
				throw Exception( "Unable to create instance of matrix type " + dist->Name(), FileAndLine );
			}
		}

		//void GetDifferenceDistributions(
		//	Histogram<Byte>& symbolDistribution,
		//	int maxK,
		//	vector<pair<int, IntegerDistribution>>& distributions
		//) {
		//	Histogram<int> h;
		//	h.GetOneMerHistogram<Byte, int>( symbolDistribution, [=]( Byte x, Byte y ) { return Difference( x, y ); } );

		//	assert_intsEqual( 0, h.data.begin()->first );
		//	assert_intsEqual( maxValue - minValue, h.data.rbegin()->first );

		//	int min = 0;
		//	int max = maxValue - minValue;

		//	vector<double> p( 1 + max - min );

		//	for ( auto kvp : h.data ) {
		//		p[(int) kvp.first - min] = kvp.second;
		//	}

		//	IntegerDistribution oneMerDistances( min, max, p );
		//	IntegerDistribution* latestDistribution = &oneMerDistances;
		//	distributions.push_back( pair<int, IntegerDistribution>( 1, oneMerDistances ) );

		//	for ( int k = 2; k <= maxK; k++ ) {
		//		IntegerDistribution nextDistribution = oneMerDistances.Add( *latestDistribution );
		//		distributions.push_back( pair<int, IntegerDistribution>( k, nextDistribution ) );
		//		latestDistribution = &oneMerDistances;
		//	}
		//}

		/**
		 *	<summary>
		 *		Fills the "top-left" corner of the supplied matrix with corresponding distance
		 *		values.
		 *	</summary>
		 *	<returns>
		 *
		 *	</returns>
		 */
		void PopulateDistanceTable( Distance lookup[128][128] ) {
			const string& symbols = alphabet->Symbols();

			for ( size_t i = 0; i < symbols.length(); i++ ) {
				for ( size_t j = 0; j < symbols.length(); j++ ) {
					lookup[i][j] = maxValue - dict[alphabet->Encode( symbols[i] ).value][alphabet->Encode( symbols[j] ).value];
				}
			}
		}

		QutBio::Alphabet* Alphabet() const {
			return alphabet;
		}

		SimilarityMatrix& SetAlphabet( QutBio::Alphabet* alphabet ) {
			this->alphabet = alphabet;
			return *this;
		}
	};

	typedef SimilarityMatrix* pSimilarityMatrix;
}
