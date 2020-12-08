#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>
#include "Array.hpp"

using namespace std;

namespace QutBio{

	class Simproj {
	public:

		/// <summary>
		/// Similarity projection using substitution distance with query 
		///	step size 1.
		/// </summary>
		/// <typeparam name="SymType">The underlying symbol type, probably byte for DNA or Amino Acid.</typeparam>
		/// <typeparam name="SymDistanceFunction">A type which defines opertor()(SymType, SymType) -> int.</typeparam>
		/// <param name="q">Query sequence.</param>
		/// <param name="r">Reference Sequence</param>
		/// <param name="k">Word length</param>
		/// <param name="symDist">Substitution matrix.</param>
		/// <param name="rowMinima">Vector which will be overwritten with row minima</param>
		/// <param name="colMinima">Vector which will be overwritten with column minima</param>
		template<typename SymType, typename SymDistanceFunction, typename Processor>
		static void ComputeKmerDistances(
			const vector<SymType>& q,
			const vector<SymType>& r,
			size_t k,
			const SymDistanceFunction& symDist,
			Processor process
		) {
			auto updateUpperTriang = [=]( size_t i, size_t j, int d ) {
				process(i,j,d);
			};

			auto updateLowerTriang = [=]( size_t i, size_t j, int d ) {
				process(j,i,d);
			};

			Triang<SymType, SymDistanceFunction, decltype( updateUpperTriang )>( q, r, 0, k, symDist, updateUpperTriang );
			Triang<SymType, SymDistanceFunction, decltype( updateLowerTriang )>( r, q, 1, k, symDist, updateLowerTriang );
		}

		/// <summary>
		/// Similarity projection using substitution distance with query 
		///	step size 1.
		/// </summary>
		/// <typeparam name="SymType">The underlying symbol type, probably byte for DNA or Amino Acid.</typeparam>
		/// <typeparam name="SymDistanceFunction">A type which defines opertor()(SymType, SymType) -> int.</typeparam>
		/// <param name="q">Query sequence.</param>
		/// <param name="r">Reference Sequence</param>
		/// <param name="k">Word length</param>
		/// <param name="symDist">Substitution matrix.</param>
		/// <param name="rowMinima">Vector which will be overwritten with row minima</param>
		/// <param name="colMinima">Vector which will be overwritten with column minima</param>
		template<typename SymType, typename SymDistanceFunction>
		static void ComputeKmerDistances(
			const vector<SymType>& q,
			const vector<SymType>& r,
			size_t k,
			const SymDistanceFunction& symDist,
			vector<int>& rowMinima,
			vector<int>& colMinima
		) {
			if ( rowMinima.size() < q.size() + 1 - k || colMinima.size() < r.size() + 1 - k ) {
				throw ArgumentException( "insufficient room in row- and/or column minimum vector.", FileAndLine );
			}

			std::fill( rowMinima.begin(), rowMinima.end(), numeric_limits<int>::max() );
			std::fill( colMinima.begin(), colMinima.end(), numeric_limits<int>::max() );

			auto updateUpperTriang = [&rowMinima, &colMinima]( size_t i, size_t j, int d ) {
				if ( d < rowMinima[i] ) rowMinima[i] = d;
				if ( d < colMinima[j] ) colMinima[j] = d;
			};

			auto updateLowerTriang = [&rowMinima, &colMinima]( size_t i, size_t j, int d ) {
				if ( d < rowMinima[j] ) rowMinima[j] = d;
				if ( d < colMinima[i] ) colMinima[i] = d;
			};

			Triang<SymType, SymDistanceFunction, decltype( updateUpperTriang )>( q, r, 0, k, symDist, updateUpperTriang );
			Triang<SymType, SymDistanceFunction, decltype( updateLowerTriang )>( r, q, 1, k, symDist, updateLowerTriang );
		}

	/// <summary>
	/// Similarity projection using substitution distance with query 
	///	step size 1.
	/// </summary>
	/// <typeparam name="SymType">The underlying symbol type, probably bye for DNA or Amino Acid.</typeparam>
	/// <param name="q">Query sequence.</param>
	/// <param name="r">REference Sequence</param>
	/// <param name="k">Word length</param>
	/// <param name="matrix">Substitution matrix.</param>
	/// <param name="rowMinima">Vector which will be overwritten with row minima</param>
	/// <param name="colMinima">Vector which will be overwritten with column minima</param>
	/// <param name="dist">Matrix which will be populated with pairwise kmer distances</param>
	template<typename SymType, typename SymDistanceFunction>
	static void ComputeKmerDistances(
		const vector<SymType>& q,
		const vector<SymType>& r,
		size_t k,
		const SymDistanceFunction& symDist,
		vector<int>& rowMinima,
		vector<int>& colMinima,
		FlatMatrix<int>& dist
	) {
		if ( rowMinima.size() < q.size() + 1 - k || colMinima.size() < r.size() + 1 - k ) {
			throw ArgumentException( "insufficient room in row- and/or column minimum vector.", FileAndLine );
		}

		std::fill( rowMinima.begin(), rowMinima.end(), numeric_limits<int>::max() );
		std::fill( colMinima.begin(), colMinima.end(), numeric_limits<int>::max() );

		auto updateUpperTriang = [&dist, &rowMinima, &colMinima]( size_t i, size_t j, int d ) {
			dist( i, j ) = d;
			if ( d < rowMinima[i] ) rowMinima[i] = d;
			if ( d < colMinima[j] ) colMinima[j] = d;
		};

		auto updateLowerTriang = [&dist, &rowMinima, &colMinima]( size_t i, size_t j, int d ) {
			dist( j, i ) = d;
			if ( d < rowMinima[j] ) rowMinima[j] = d;
			if ( d < colMinima[i] ) colMinima[i] = d;
		};

		Triang<SymType, SymDistanceFunction, decltype( updateUpperTriang )>( q, r, 0, k, symDist, updateUpperTriang );
		Triang<SymType, SymDistanceFunction, decltype( updateLowerTriang )>( r, q, 1, k, symDist, updateLowerTriang );
	}

	template<typename SymType, typename SymDistanceFunction, typename Update>
	static void Triang(
		const vector<SymType>& q,
		const vector<SymType>& r,
		size_t offset, /* should be zero for upper triangle, or 1 for lower triangle */
		size_t k,
		const SymDistanceFunction& symDist,
		Update update
	) {
		if ( k > q.size() ) return;
		if ( k > r.size() ) return;

		const size_t m = q.size() + 1 - k;
		const size_t n = r.size() + 1 - k;

		auto buffer = (int*) alloca( k * sizeof( int ) );
		int d = 0;

		for ( size_t c = offset; c < n; c++ ) {
			int d = 0;
			size_t t = 0;

			for ( ; t < k; t++ ) {
				buffer[t] = symDist( q[t], r[c + t] );
				d += buffer[t];
			}

			size_t i = 0, j = c;
			update( i, j, d );
			i++;
			j++;

			for ( ; i < m && j < n; i++, j++, t++ ) {
				size_t nextPos = t % k;
				d -= buffer[nextPos];
				auto dt = symDist( q[t], r[c + t] );
				buffer[nextPos] = dt;
				d += dt;
				update( i, j, d );
			}
		}
	}

	template<typename T>
	static double BestOfBest( const vector<T>& rowMinima, size_t rowCount, const vector<T>& colMinima, size_t colCount ) {
		double minimum = rowMinima[0];

		for ( size_t i = 1; i < rowCount; i++ ) {
			double d = rowMinima[i];

			if ( d < minimum ) {
				minimum = d;
			}
		}

		for ( size_t i = 0; i < colCount; i++ ) {
			double d = colMinima[i];

			if ( d < minimum ) {
				minimum = d;
			}
		}

		return minimum;
	}

	template<typename T>
	static double HausdorffAverageAverage( const vector<T>& rowMinima, size_t rowCount, const vector<T>& colMinima, size_t colCount ) {
		double rowTot = rowMinima[0];

		for ( size_t i = 1; i < rowCount; i++ ) {
			rowTot += rowMinima[i];
		}

		double colTot = colMinima[0];

		for ( size_t i = 1; i < colCount; i++ ) {
			colTot += colMinima[i];
		}

		return ( rowTot / rowCount + colTot / colCount ) / 2;
	}

	/// <summary>
	/// Gets the modified Hausdorff distance defined by Jain & Dubuisson from a pair of row- and column minimum vectors.
	/// </summary>
	/// <typeparam name="T">The element type in row- and column minima.</typeparam> 
	/// <param name="rowMinima">Vector containing the row minima in cells [0,rowCount).</param>
	/// <param name="rowCount">Number of meaningful entries in rowMinima.</param>
	/// <param name="colMinima">Vector containing the column minima in cells [0,colCount).</param>
	/// <param name="colCount">Number of meaningful entries in colMinima.</param>
	/// <returns>Returns max( avg( rowMinima[0,rowCount) ), avg( colMinima[0,colCount) ) )</returns>

	template<typename T>
	static double HausdorffAverage( const vector<T>& rowMinima, size_t rowCount, const vector<T>& colMinima, size_t colCount ) {
		double rowTot = rowMinima[0];

		for ( size_t i = 1; i < rowCount; i++ ) {
			rowTot += rowMinima[i];
		}

		double colTot = colMinima[0];

		for ( size_t i = 1; i < colCount; i++ ) {
			colTot += colMinima[i];
		}

		return std::max( rowTot / rowCount, colTot / colCount );
	}

	/// <summary>
	/// Gets the "classic" Hausdorff distance from a pair of row- and column minimum vectors.
	/// </summary>
	/// <param name="rowMinima">Vector containing the row minima in cells [0,rowCount).</param>
	/// <param name="rowCount">Number of meaningful entries in rowMinima.</param>
	/// <param name="colMinima">Vector containing the column minima in cells [0,colCount).</param>
	/// <param name="colCount">Number of meaningful entries in colMinima.</param>
	/// <returns>Returns max( max( rowMinima[0,rowCount) ), max( colMinima[0,colCount) ) )</returns>
	
	template<typename T>
	static double Hausdorff( const vector<T>& rowMinima, size_t rowCount, const vector<T>& colMinima, size_t colCount ) {
		double rowMax = rowMinima[0];

		for ( size_t i = 1; i < rowCount; i++ ) {
			if ( rowMinima[i] > rowMax ) rowMax = rowMinima[i];
		}

		double colMax = colMinima[0];

		for ( size_t i = 1; i < colCount; i++ ) {
			if ( colMinima[i] > colMax ) colMax = colMinima[i];
		}

		return std::max( rowMax, colMax );
	}

};

}