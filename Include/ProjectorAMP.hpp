#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "Alphabet.hpp"
#include "SimilarityMatrix.hpp"
#include "FastaSequence.hpp"
#include "FreeList.hpp"
#include "Util.hpp"
#include "SequenceDistanceFunction.hpp"

#include <string>
#include <ctime>
#include <cfloat>

#include <amp.h>
#include <iostream> 

#ifdef MARKERS
#include <cvmarkersobj.h>
using namespace concurrency::diagnostic;
marker_series g_markerSeries(L"KmerRankCpp");
#endif

using namespace concurrency;

#undef max
#undef min

namespace QutBio {

	class CachedSequenceArrayView {
	public:
		const char * charData;
		vector<int> cpuData;
		array_view<const int, 1> gpuData;

		CachedSequenceArrayView(const char * charData, int len)
			: cpuData(len), gpuData(len, cpuData) {

			for ( int i = 0; i < len; i++ ) {
				cpuData[i] = charData[i];
			}
		}
	};


	class ProjectorAMP : SequenceDistanceFunction {
	protected:

		accelerator accel;
		accelerator_view view;
		array_view<const int, 2> &cd;
		array<int, 2> dist;
		array<int, 2> transposeDist;
		array_view<int, 1> & rowTotal;
		array_view<int, 1> &colTotal;

	public:
		ProjectorAMP(
			SimilarityMatrix * matrix
			, uint kmerLength
			, array_view<const int, 2> &cd
			, int m
			, int n
			, array_view<int, 1> &rowTotal
			, array_view<int, 1> &colTotal
			) :
			SequenceDistanceFunction(matrix, kmerLength)
			, cd(cd)
			, accel(accelerator::default_accelerator)
			, view(accel.default_view)
			, dist(m, n, view)
			, transposeDist(m, n, view)
			, rowTotal(rowTotal)
			, colTotal(colTotal) {}

		virtual ~ProjectorAMP() {}

		double ComputeDistance(
			EncodedFastaSequence * querySeq,
			EncodedFastaSequence * subjectSeq
			) {
			throw NotImplementedException(FileAndLine);
		}

		double ComputeDistance(
			FastaSequence * querySeq,
			EncodedFastaSequence * subjectSeq
			, array_view<int, 1> &rowTotal
			, array_view<int, 1> &colTotal
			, int subjectIdx
			) {
			static unordered_map<string, CachedSequenceArrayView *> cache;

			string &queryStr = querySeq->Sequence();
			const char * queryChars = queryStr.c_str();
			const size_t queryKmerCount = queryStr.size() - kmerLength + 1;

			string & subjectStr = subjectSeq->Sequence();
			const char * subjectChars = subjectStr.c_str();
			const size_t subjectKmerCount = subjectStr.size() - kmerLength + 1;

			const int m = (int) queryKmerCount;
			const int n = (int) subjectKmerCount;

			CachedSequenceArrayView * queryCache;
			CachedSequenceArrayView * subjectCache;

			auto q = cache.find(queryStr);

			if ( q == cache.end() ) {
				queryCache = new CachedSequenceArrayView(queryChars, (int) querySeq->Sequence().size());
				cache[queryStr] = queryCache;
			}
			else {
				queryCache = q->second;
			}

			q = cache.find(subjectStr);

			if ( q == cache.end() ) {
				subjectCache = new CachedSequenceArrayView(subjectChars, (int) subjectSeq->Sequence().size());
				cache[subjectStr] = subjectCache;
			}
			else {
				subjectCache = q->second;
			}

			ComputeDistanceProfilesAMP(queryCache->gpuData, subjectCache->gpuData, cd, rowTotal, colTotal, subjectIdx,
				pair<array<int, 2> &, array<int, 2> &>(dist, transposeDist)
				);

			return 0;

		}

		/// <summary>
		///	Given strings s and t, having length m and n respectively, and a similarity matrix 
		///	reformatted as a 128 by 128 array of integers, populates a vector contain a kmer mutual
		///	distance table.
		/// </summary>
		///	<param name="queryBytes"></param>

	public:
		// Constructs a character mutual distance table from a similarity matrix. 
		static void GetCharDistMatrix(
			const SimilarityMatrix * m,
			int charDist[]
			) {
			for ( int i = 0; i < 128; i++ ) {
				for ( int j = 0; j < 128; j++ ) {
					int dist = m->Difference((char) i, (char) j);

					charDist[i * 128 + j] = dist;
				}
			}
		}

	private:
		void ComputeDistanceProfilesAMP(
			array_view<const int, 1> &query,
			array_view<const int, 1> &subject,
			array_view<const int, 2> &cd,
			array_view<int, 1> &rowTotal,
			array_view<int, 1> &colTotal,
			int subjectIdx,
			pair<array<int, 2>&, array<int, 2>&> dist
			) {
			const int k = kmerLength;
			int m = query.extent[0] - k + 1;
			int n = subject.extent[0] - k + 1;

			static bool deja = false;

			if ( !deja ) {
				// std::wcout << "Using accelerator " << view.accelerator.description << endl;
				deja = true;
			}

			// g_markerSeries.write_flag(diagnostic:: normal_importance, L" Create arrays");

			//array<int, 2> dist.first(m,n, view );
			//array<int, 2> dist.second(m,n, view );

			// g_markerSeries.write_flag(diagnostic:: normal_importance, L"Begin distance calc");

			parallel_for_each(view, concurrency::extent<2>(m, n),
				[= /*, &dist.first, &dist.second */](index<2> idx) restrict(amp)
			{
				int row = idx[0];
				int col = idx[1];
				int sum = 0;
				for ( int i = 0; i < k; i++ ) {
					int x = query[row + i];
					int y = subject[col + i];
					sum += cd[x][y];
				}
				dist.first[idx] = sum;
				dist.second[idx] = sum;
			}
			);
			// g_markerSeries.write_flag(diagnostic:: normal_importance, L"End distance calc");

#if 0
			ofstream f("a:/Temp/debug_output.txt");
			f << "Distances" << endl;
			dumpDistanceArray(f, m, n, distance.first, queryChars, subjectChars);
#endif // 1

			// Get rowMinima in column 0 of dist.first
			// Implementation based on AMP book, but modified to eliminate
			// a nasty bug.
			{
				int N = n;
				const int minItems = 1;

				while ( N > minItems ) {
					// N is Number of items remaining to be processed in each iteration.
					// s is the number of items that will be folded back from high indices to lower indices.
					// max is the highest relevant index to be processed.
					// excess is the number of items that will be remaining to process at the end of this iteration.

					int s = N >> 1;
					int max = N - 1;  // 
					int excess = N - s;
					concurrency::extent<2> grid(m, s);

					parallel_for_each(view, grid, [= /*, &dist.first*/](index<2> idx) restrict(amp)
					{
						int row = idx[0];
						int srcCol = max - idx[1];
						int destCol = srcCol - excess;

						auto d1 = dist.first[row][destCol];
						auto d2 = dist.first[row][srcCol];

						if ( d2 < d1 ) {
							dist.first[row][destCol] = d2;
						}
					});

					N = excess;

#if 0
					f << "After stride = " << s << ":" << endl;
					dumpDistanceArray(f, m, n, distance.first, queryChars, subjectChars);
#endif // 0
				}

				//// Clean up the last few columns
				//concurrency::extent<1> grid2(m);

				//parallel_for_each(view, grid2, [=, &dist.first](index<1> idx) restrict(amp)
				//{
				//	int row = idx[0];

				//	for ( int col = 1; col < N; col++ ) {
				//		auto d1 = dist.first[row][0];
				//		auto d2 = dist.first[row][col];

				//		if ( d2 < d1 ) {
				//			dist.first[row][col] = d2;
				//		}
				//	}
				//});
			}

#if 0
			vector<int> dbRowMin(m);
			copy(distance.first.section(0, 0, m,1), dbRowMin.begin());

			f << endl;

			for ( int row = 0; row < m; row++ ) {
				f << queryChars[row] << "\t" << dbRowMin[row] << endl;
			}

			f << endl;
#endif // 0


			// Get rowMinima total in dist.first[0][0]
			for ( int N = m; N > 1; ) {
				int s = N >> 1;
				int max = N - 1;  // 
				int excess = N - s;
				concurrency::extent<1> grid(s);

				parallel_for_each(view, grid, [= /*, &dist.first*/](index<1> idx) restrict(amp)
				{
					int srcRow = max - idx[0];
					int destRow = srcRow - excess;
					dist.first[destRow][0] += dist.first[srcRow][0];
				});

				N = excess;
			}

#if 0
			vector<int> dbRowMinTot(1);
			copy(distance.first.section(0, 0, 1, 1), dbRowMinTot.begin());

			f << endl;

			f << "tot" << "\t" << dbRowMinTot[0] << endl;

			f << endl;

			f.flush();
			f.close();
			abort();
#endif // 0


			// Get colMinima in dist.second[0][*]
			for ( int N = m; N > 1; ) {
				int s = N >> 1;
				int max = N - 1;
				int excess = N - s;
				concurrency::extent<2> grid(s, n);

				parallel_for_each(view, grid, [=/*, &dist.second*/](index<2> idx) restrict(amp)
				{
					int srcRow = max - idx[0];
					int destRow = srcRow - excess;
					int col = idx[1];

					auto d1 = dist.second[destRow][col];
					auto d2 = dist.second[srcRow][col];

					if ( d2 < d1 ) {
						dist.second[destRow][col] = d2;
					}
				});

				N = excess;

#if 0
				f << "After stride = " << s << ":" << endl;
				dumpDistanceArray(f, m, n, distance.first, queryChars, subjectChars);
#endif // 0
			}

			// Get colMinima total in dist.second[0][0]
			for ( int N = n; N > 1; ) {
				int s = N >> 1;
				int max = N - 1;
				int excess = N - s;
				concurrency::extent<1> grid(s);

				parallel_for_each(view, grid, [=/*, &dist.second*/](index<1> idx) restrict(amp)
				{
					int srcCol = max - idx[0];
					int destCol = srcCol - excess;
					dist.second[0][destCol] += dist.second[0][srcCol];
				});

				N = excess;
			}

			// Hacky append results to row and column total vectors.
			{
				concurrency::extent<1> grid(1);

				parallel_for_each(view, grid,
					[=/*, &dist.first, &dist.second*/](index<1> idx) restrict(amp)
				{
					rowTotal[subjectIdx] = dist.first[0][0];
					colTotal[subjectIdx] = dist.second[0][0];
				});
			}
		}

		static void dumpDistanceArray(ostream & f, int m, int n, array<int, 2> & dist1, const char * queryChars, const char * subjectChars) {
			vector<int> dbDist(m * n);
			array_view<int, 2> dbDistCpu(m, n, dbDist);
			copy(dist1, dbDist.begin());

			f << "\t";

			for ( int col = 0; col < n; col++ ) {
				f << "\t" << subjectChars[col];
			}

			f << endl;

			for ( int row = 0; row < m; row++ ) {
				f << queryChars[row] << "\t";

				for ( int col = 0; col < n; col++ ) {
					f << "\t" << dbDistCpu[row][col];
				}

				f << endl;
			}
		}

	};
}
