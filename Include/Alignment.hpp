#pragma once
#include <functional>

using namespace std;

#include "Array.hpp"
#include "Types.hpp"

namespace QutBio {
	class Alignment {
	public:
		enum Transition { Zero, Horizontal, Vertical, Diagonal };
		enum Mode { Needle, Water };

	protected:
		size_t M;
		size_t N;
		FlatMatrix<int> H;
		FlatMatrix<int> P;
		FlatMatrix<int> Q;
		FlatMatrix<Transition> T;
		size_t IMax;
		size_t JMax;

	public:

		/// <summary>
		/// The length of the "query" sequence.
		/// </summary>
		size_t GetM() const { return this->M; }

		/// <summary>
		/// The length of the "reference" sequence.
		/// </summary>
		size_t GetN() const { return this->N; }

		/// <summary>
		/// The cumulative path score matrix. This is a (m+1)\times(n+1) matrix.
		/// H[0,j] = H[i,0] = 0 for 1 <= i <= m and 1 <= j <= n. 
		/// H[i,j] contains the cumulative score reached by travelling from 
		/// (a[0],b[0]) to (a[i-1],b[j-1]) along an optimal alignment.
		/// </summary>
		const FlatMatrix<int>& GetH() const { return H; }

		/// <summary>
		/// Auxiliary storage for vertical gap penalty.
		/// </summary>
		const FlatMatrix<int>& GetP() const { return this->P; }

		/// <summary>
		/// Auxiliary storage for horizontal gap penalty.
		/// </summary>
		const FlatMatrix<int>& GetQ() const { return this->Q; }

		/// <summary>
		/// Record of transitions. At each step four possibilities are present.
		/// Zero - Score resets to zero (Smith-Waterman only)
		/// Horizontal - Vertical gap increment; horizontal move in matrix.
		/// Vertical - Horizontal gap increment; vertical move in matrix.
		/// Diagonal - Accept current match; diagonal move in matrix.
		/// </summary>
		const FlatMatrix<Transition>& GetT() const { return this->T; }

		/// <summary>
		/// The position in the "query" sequence at which max(H[i,j]) occurs. 
		/// </summary>
		size_t GetIMax() const { return this->IMax; }

		/// <summary>
		/// The position in the "reference" sequence at which max(H[i,j]) occurs. 
		/// </summary>
		size_t GetJMax() const { return this->JMax; }


		/// <summary>
		/// 
		/// </summary>
		/// <param name="a">An encoded sequence to be compared (the "query"). Each symbol is a sing</param>
		/// <param name="b">An encoded sequence to be compared (the "reference").</param>
		/// <param name="sim">A callable object which returns the pairwise similarity between two symbols.</param>
		/// <param name="u">u in the affine gap penalty: g(k) = uk + v.</param>
		/// <param name="v">v in the affine gap penalty: g(k) = uk + v.</param>
		template <typename SymbolType, typename SimilarityFunc>
		void Align( 
			const vector<SymbolType> & a,
			const vector<SymbolType> & b, 
			SimilarityFunc sim, 
			int u, 
			int v, 
			Mode alignmentMode 
		) {
			M = a.size();
			N = b.size();

			H.resize(M + 1, N + 1);
			P.resize(M + 1, N + 1);
			Q.resize(M + 1, N + 1);
			T.resize(M + 1, N + 1);
			IMax = JMax = 0;

			for ( size_t i = 1; i <= M; i++ ) {
				for ( size_t j = 1; j <= N; j++ ) {
					P(i, j) = std::max( H(i - 1, j) - v, P(i - 1, j) ) - u;
					Q(i, j) = std::max( H(i, j - 1) - v, Q(i, j - 1) ) - u;

					H(i, j) = H(i - 1, j - 1) + sim( a[i - 1], b[j - 1] );
					T(i, j) = Transition::Diagonal;

					if ( P(i, j) > H(i, j) ) {
						H(i, j) = P(i, j);
						T(i, j) = Transition::Vertical;
					}

					if ( Q(i, j) > H(i, j) ) {
						H(i, j) = Q(i, j);
						T(i, j) = Transition::Horizontal;
					}

					if ( alignmentMode == Mode::Water ) {
						if ( 0 > H(i, j) ) {
							H(i, j) = 0;
							T(i, j) = Transition::Zero;
						}

						if ( H(i, j) > H(IMax, JMax) ) {
							IMax = i;
							JMax = j;
						}
					}
					else {
						if ( ( i == M || j == N ) && ( H(i, j) > H(IMax, JMax) ) ) {
							IMax = i;
							JMax = j;
						}
					}
				}
			}
		}

		auto MaxScore() const {
			return H.at(IMax, JMax);
		}

		struct Coord {
			size_t I;
			size_t J;
			Transition T;
			friend ostream & operator<<( ostream & w, const Coord & coord ) {
				w << coord.I << "\t" << coord.J << "\t" << coord.T;
				return w;
			}
		};

		vector<Coord> TraceFrom( size_t i, size_t j ) {
			vector<Coord> path;

			while ( T(i, j) > Transition::Zero && i > 0 && j > 0 ) {
				Transition t = T(i, j);
				Coord c{i, j, t};
				path.push_back(c);

				switch ( t ) {
				case Transition::Vertical: i--; break;
				case Transition::Horizontal: j--; break;
				default: i--; j--; break;
				}
			}

			std::reverse(path.begin(), path.end());
			return path;
		}
	};
}

