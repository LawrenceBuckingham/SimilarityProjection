#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "Exception.hpp"
#include "Subrange.hpp"
#include "Types.hpp"

namespace QutBio {

	using namespace std;

	class SubstitutionMatrix {
	private:
		vector<string> metadata;
		string alphabet;
		byte inverse[256];
		size_t size;
		int min, max;

		static const size_t MAT_SIZE = 256;
		static const byte MISSING = MAT_SIZE - 1;

		int similarity[MAT_SIZE][MAT_SIZE];
		int distance[MAT_SIZE][MAT_SIZE];

		FlatMatrix<int> digramDist;

	public:
		static Symbol Missing() {
			Symbol missing{ (byte) MISSING };
			return missing;
		}

		/// <summary>
		/// Returns the numeric symbol corresponding to a given char value.
		///	If the char value is not in the domain of this substitution matrix, 
		///	the result will be SubstitutionMstrix::Missing(). 
		/// </summary>
		/// <param name="ch"></param>
		/// <returns></returns>
		Symbol Encode( char ch ) const {
			Symbol res{ inverse[(byte) ch] };
			return res;
		}

		bool IsDefined( char ch ) const {
			return inverse[byte( ch )] != MISSING;
		}

		vector<string>& Metadata() { return metadata; }

		const string& Alphabet() const { return alphabet; }

		vector<Symbol> Symbols() const {
			vector<Symbol> symbols;
			for ( byte i = 0; i < (byte)(size); i++ ) {
				symbols.push_back(Symbol{i});
			}
			return symbols;
		}

		int Similarity( Symbol x, Symbol y ) const {
			return similarity[x.value][y.value];
		}

		size_t Size() const { return size; }

		int operator()( Symbol x, Symbol y ) const { return distance[x.value][y.value]; }

		int operator()( Digram x, Digram y ) const { return digramDist.at( x, y ); }

		int operator()( const Symbol* x, const Symbol* y, size_t k ) const {
			int d = 0;

			for ( size_t i = 0; i < k; i++ ) {
				d += distance[x[i].value][y[i].value];
			}

			return d;
		}

		int operator()( const Digram* x, const Digram* y, size_t k ) const {
			int d = 0;

			for ( size_t i = 0; i < k; i += 2 ) {
				d += digramDist.at( x[i], y[i] );
			}

			return d;
		}

		int Min() const { return min; }

		int Max() const { return max; }

		SubstitutionMatrix( const std::string& substitutionMatrix ) {
			istringstream r( substitutionMatrix );
			Parse( r );
		}

		SubstitutionMatrix( std::istream& stream ) {
			Parse( stream );
		}

		void Parse( std::istream& r ) {
			memset( inverse, 255, sizeof( inverse ) );

			bool waitingForHeadings = true;
			string line;
			int row = 0;

			while ( std::getline( r, line ) ) {
				if ( line[0] == '#' ) {
					metadata.push_back( line );
					continue;
				}

				auto fields = String::Split( line, ", \t" );

				if ( waitingForHeadings ) {
					waitingForHeadings = false;

					for ( auto& f : fields ) {
						alphabet.push_back( f[0] );
					}

					size = alphabet.size();

					for ( size_t i = 0; i < size; i++ ) {
						inverse[tolower( alphabet[i] )] = i;
						inverse[toupper( alphabet[i] )] = i;
					}

					for ( size_t i = 0; i < MAT_SIZE; i++ ) {
						for ( size_t j = 0; j < MAT_SIZE; j++ ) {
							similarity[i][j] = numeric_limits<int>::min();
							distance[i][j] = numeric_limits<int>::max();
						}
					}
				}
				else {
					int offset = 0;
					int idx = row;

					if ( fields.size() < size ) {
						break;
					}
					else if ( fields.size() > size ) {
						try {
							offset = 1;
							idx = inverse[ fields[0][0] ];
						}
						catch ( out_of_range ) {
							idx = -1;
						}
					}

					if ( idx < 0 || idx >= (int) size ) {
						cerr << "idx out of bounds for symbol '" << fields[0] << "': " << idx << " -- input line skipped\n";
					}
					else {
						for ( size_t col = 0; col < size; col++ ) {
							int sim = Int::Parse( fields[col + offset] );
							similarity[idx][col] = sim;
							if ( sim < min ) min = sim;
							if ( sim > max ) max = sim;
						}
					}

					row++;
				}
			}

			for ( size_t i = 0; i < size; i++ ) {
				for ( size_t j = 0; j < size; j++ ) {
					distance[i][j] = max - similarity[i][j];
				}
			}

			this->digramDist.resize( size * size, size * size );

			for ( size_t i = 0; i < size; i++ ) {
				for ( size_t j = 0; j < size; j++ ) {
					for ( size_t k = 0; k < size; k++ ) {
						for ( size_t l = 0; l < size; l++ ) {
							digramDist( i + j * size, k + l * size ) = distance[i][k] + distance[j][l];
						}
					}
				}
			}
		}

		friend ostream& operator << ( ostream& out, const SubstitutionMatrix& mat ) {
			for ( auto s : mat.metadata ) out << s << "\n";

			bool deja = false;

			for ( auto c : mat.alphabet ) {
				out << ( deja ? "\t" : "" ) << c;
				deja = true;
			}

			out << "\n";

			for ( auto c : mat.alphabet ) {
				auto x = mat.Encode( c );
				out << ( c );

				deja = false;

				for ( auto d : mat.alphabet ) {
					auto y = mat.Encode( d );
					out << ( deja ? "\t" : "" ) << mat.similarity[x.value][y.value];
					deja = true;
				}

				out << "\n";
			}

			return out;
		}

		/// <summary>
		/// Maps symbols in the alphabet supported by this matrix into the half-open 
		/// subrange [0,size), wrapping the results in Ordinal<byte> for type safety.
		/// </summary>
		/// <param name="s"></param>
		/// <returns></returns>
		template<typename iterator>
		void Encode( const iterator& begin, const iterator& end, vector<Symbol>& res ) const {
			res.resize( end - begin );
			int i = 0;

			for ( auto p = begin; p != end; p++ ) {
				res[i++] = Encode( (char) *p );
			}
		}

		char Decode(Symbol symbol) const {
			return alphabet[symbol.value];
		}

		string Decode( const vector<Symbol>& seq ) const {
			string s( seq.size(), ' ' );
			int i = 0;

			for ( auto x : seq ) {
				s[i++] = alphabet[x.value];
			}

			return s;
		}

		int MinSim() const { return min; }

		int MaxSim() const { return max; }

		int MinDist() const { return 0; }

		int MaxDist() const { return max - min; }

		static const SubstitutionMatrix& Blosum62() {
			static SubstitutionMatrix instance{ _blosum62() };
			return instance;
		}

		static const SubstitutionMatrix& Fitch() {
			static SubstitutionMatrix instance{ _fitch() };
			return instance;
		}

		static string _blosum62() {
#pragma region _blosum62
			static string s =
				"#  Matrix made by matblas from blosum62.iij"  "\n"
				"#  * column uses minimum score"  "\n"
				"#  BLOSUM Clustered Scoring Matrix in 1 / 2 Bit Units"  "\n"
				"#  Blocks Database = / data / blocks_5.0 / blocks.dat"  "\n"
				"#  Cluster Percentage : >= 62"  "\n"
				"#  Entropy = 0.6979, Expected = -0.5209"  "\n"
				"A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X *"  "\n"
				"A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4"  "\n"
				"R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4"  "\n"
				"N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4"  "\n"
				"D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4"  "\n"
				"C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4"  "\n"
				"Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4"  "\n"
				"E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4"  "\n"
				"G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4"  "\n"
				"H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4"  "\n"
				"I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4"  "\n"
				"L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4"  "\n"
				"K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4"  "\n"
				"M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4"  "\n"
				"F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4"  "\n"
				"P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4"  "\n"
				"S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4"  "\n"
				"T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4"  "\n"
				"W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4"  "\n"
				"Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4"  "\n"
				"V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4"  "\n"
				"B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4"  "\n"
				"Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4"  "\n"
				"X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4"  "\n"
				"* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1"  "\n"
				"";
			return s;
#pragma endregion
		}

		static string _fitch() {
#pragma region _fitch
			static string s =
				"   D  C  T  F  E  H  K  A  M  N  Y  P  Q  R  S  W  L  V  I  G  *  B  Z  X" "\n"
				"D -0 -2 -2 -2 -1 -1 -2 -1 -3 -1 -1 -2 -2 -2 -2 -3 -2 -1 -2 -1 -3 -3 -2 -1" "\n"
				"C -2 -0 -2 -1 -3 -2 -3 -2 -3 -2 -1 -2 -2 -1 -1 -1 -2 -1 -2 -1 -3 -1 -1 -1" "\n"
				"T -2 -2 -0 -2 -2 -2 -1 -1 -1 -1 -2 -1 -2 -1 -1 -2 -2 -2 -1 -2 -3 -1 -1 -2" "\n"
				"F -2 -1 -2 -0 -3 -2 -3 -2 -2 -2 -1 -2 -2 -2 -1 -2 -1 -1 -1 -2 -3 -1 -1 -1" "\n"
				"E -1 -3 -2 -3 -0 -2 -1 -1 -2 -2 -2 -2 -1 -2 -2 -2 -2 -1 -3 -1 -3 -2 -3 -1" "\n"
				"H -1 -2 -2 -2 -2 -0 -2 -2 -3 -1 -1 -1 -1 -1 -2 -3 -1 -2 -2 -2 -3 -2 -2 -1" "\n"
				"K -2 -3 -1 -3 -1 -2 -0 -2 -1 -1 -2 -2 -1 -1 -2 -2 -2 -2 -2 -2 -3 -2 -2 -1" "\n"
				"A -1 -2 -1 -2 -1 -2 -2 -0 -2 -2 -2 -1 -2 -2 -1 -2 -2 -1 -2 -1 -3 -1 -1 -2" "\n"
				"M -3 -3 -1 -2 -2 -3 -1 -2 -0 -2 -3 -2 -2 -1 -2 -2 -1 -1 -1 -2 -3 -1 -1 -1" "\n"
				"N -1 -2 -1 -2 -2 -1 -1 -2 -2 -0 -1 -2 -2 -2 -1 -3 -2 -2 -1 -2 -3 -2 -2 -1" "\n"
				"Y -1 -1 -2 -1 -2 -1 -2 -2 -3 -1 -0 -2 -1 -2 -1 -2 -2 -2 -2 -2 -3 -1 -1 -1" "\n"
				"P -2 -2 -1 -2 -2 -1 -2 -1 -2 -2 -2 -0 -1 -1 -1 -2 -1 -2 -2 -2 -3 -1 -1 -1" "\n"
				"Q -2 -2 -2 -2 -1 -1 -1 -2 -2 -2 -1 -1 -0 -1 -1 -1 -1 -2 -3 -2 -3 -2 -2 -1" "\n"
				"R -2 -1 -1 -2 -2 -1 -1 -2 -1 -2 -2 -1 -1 -0 -1 -1 -1 -2 -2 -1 -3 -1 -2 -1" "\n"
				"S -2 -1 -1 -1 -2 -2 -2 -1 -2 -1 -1 -1 -1 -1 -0 -1 -1 -2 -1 -1 -3 -2 -2 -2" "\n"
				"W -3 -1 -2 -2 -2 -3 -2 -2 -2 -3 -2 -2 -1 -1 -1 -0 -1 -2 -3 -1 -3 -1 -1 -1" "\n"
				"L -2 -2 -2 -1 -2 -1 -2 -2 -1 -2 -2 -1 -1 -1 -1 -1 -0 -1 -1 -1 -3 -1 -1 -1" "\n"
				"V -1 -1 -2 -1 -1 -2 -2 -1 -1 -2 -2 -2 -2 -2 -2 -2 -1 -0 -1 -1 -3 -1 -1 -1" "\n"
				"I -2 -2 -1 -1 -3 -2 -2 -2 -1 -1 -2 -2 -3 -2 -1 -3 -1 -1 -0 -2 -3 -1 -1 -1" "\n"
				"G -1 -1 -2 -2 -1 -2 -2 -1 -2 -2 -2 -2 -2 -1 -1 -1 -1 -1 -2 -0 -3 -1 -1 -1" "\n"
				"* -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -3 -1 -3 -2 -1" "\n"
				"B -3 -1 -1 -1 -2 -2 -2 -1 -1 -2 -1 -1 -2 -1 -2 -1 -1 -1 -1 -1 -3 -3 -2 -1" "\n"
				"Z -2 -1 -1 -1 -3 -2 -2 -1 -1 -2 -1 -1 -2 -2 -2 -1 -1 -1 -1 -1 -2 -2 -3 -1" "\n"
				"X -1 -1 -2 -1 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1" "\n"
				"";
			return s;
#pragma endregion
		}
	};

}
