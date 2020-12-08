#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "Alphabet.hpp"
#include "SimilarityMatrix.hpp"
#include "FastaSequence.hpp"
#include "FragmentAggregationMode.hpp"
#include "FreeList.hpp"
#include "Util.hpp"
#include "Projector.hpp"

#include <string>
#include <ctime>
#include <cfloat>

namespace QutBio {

	class ProjectorSlice : public Projector {
	protected:
		uint interval;
		const bool debug = false;


	public:
		ProjectorSlice(
			SimilarityMatrix * matrix,
			uint kmerLength,
			const FragmentAggregationMode * fragMode,
			uint interval
		) :
			Projector(matrix, kmerLength, fragMode), interval(interval) {}

		virtual ~ProjectorSlice() {}

		void SetDebug(bool debug) {
			// this->debug = debug; 
		}
		bool Debug() { return this->debug; }

		/// <summary>
		///	Given strings s and t, having length m and n respectively, and a similarity matrix 
		///	reformatted as a 128 by 128 array of integers, populates a vector contain a kmer mutual
		///	distance table.
		/// </summary>
		///	<param name="queryBytes"></param>

	protected:
		virtual void ComputeDistanceMatrix(
			FastaSequence * querySeq,
			FastaSequence * subjectSeq,
			const Symbol* queryChars,
			const Symbol* subjectChars,
			vector<Distance> & rowMinima,
			vector<Distance> & colMinima,
			size_t queryKmerCount,
			size_t subjectKmerCount
		) {
			const size_t m = queryKmerCount;
			const size_t n = subjectKmerCount;

			if (m < interval || n < interval) {
				// For very short sequences, fall back on something more efficient.
				Projector::ComputeDistanceMatrix(
					// *querySeq, *subjectSeq, 
					queryChars, subjectChars, rowMinima, colMinima, queryKmerCount, subjectKmerCount);
				return;
			}

			Util::Fill(rowMinima, numeric_limits<Distance>::max());
			Util::Fill(colMinima, numeric_limits<Distance>::max());

			vector<bool> diagColsProcessed(m + n, false);
			size_t interval = this->interval;

			if (debug) cerr << "m = " << m << ", n = " << n << endl;

			if (m <= n) {
				if (interval > n / 2) {
					interval = n / 2;
				}

				if (debug) cerr << "interval = " << interval << endl;

				for (size_t c = interval; c <= n; c += interval) {
					DoSliceHighLow(queryChars, subjectChars, rowMinima, colMinima, diagColsProcessed, m, n, 0, c);
					DoSliceLowHigh(queryChars, subjectChars, rowMinima, colMinima, diagColsProcessed, m, n, m - 1, c);
				}
			}
			else {
				if (interval > m / 2) {
					interval = m / 2;
				}

				for (size_t r = interval; r <= m; r += interval) {
					DoSliceLowHigh(queryChars, subjectChars, rowMinima, colMinima, diagColsProcessed, m, n, r, 0);
					DoSliceHighLow(queryChars, subjectChars, rowMinima, colMinima, diagColsProcessed, m, n, r, n - 1);
				}
			}

			if (debug) {
				cerr << "rowMinima:" << endl;

				for (size_t i = 0; i < m; i++) {
					cerr << i << "\t" << rowMinima[i] << endl;
				}

				cerr << "colMinima:" << endl;

				for (size_t i = 0; i < n; i++) {
					cerr << i << "\t" << colMinima[i] << endl;
				}
			}
		}

		//	Cuts diagonally across the distance matrix in the direction of decreasing 
		//	query index and increasing subject index.
		void DoSliceLowHigh(
			const Symbol* queryChars,
			const Symbol* subjectChars,
			vector<Distance> & rowMinima,
			vector<Distance> & colMinima,
			vector<bool> & diagColsProcessed,
			const size_t m,
			const size_t n,
			const size_t r,
			const size_t c
		) {
			Distance minDist = MAX_DIST;
			int minRow = 0;
			int minCol = 0;

			for (int col = (int)c, row = (int)r; (size_t)col < n && row >= 0; col++, row--) {
				Distance distance = GetKmerDistance(queryChars + row, subjectChars + col, kmerLength);

				if (distance < minDist) {
					minDist = distance;
					minRow = row;
					minCol = col;
				}

				if (distance < rowMinima[row]) rowMinima[row] = distance;
				if (distance < colMinima[col]) colMinima[col] = distance;

				if (row == 0) continue;

				distance = GetKmerDistance(queryChars + row - 1, subjectChars + col, kmerLength);

				if (distance < minDist) {
					minDist = distance;
					minRow = row - 1;
					minCol = col;
				}

				if (distance < rowMinima[row - 1]) rowMinima[row - 1] = distance;
				if (distance < colMinima[col]) colMinima[col] = distance;
			}

			if (diagColsProcessed[m + minCol - minRow]) return;

			diagColsProcessed[m + minCol - minRow] = true;

			if ((size_t)minRow < m - 1 && (size_t)minCol < n - 1) {
				RunDownDiagonal(queryChars, subjectChars, rowMinima, colMinima, m, n, minRow + 1, minCol + 1);
			}

			if (minRow > 0 && minCol > 0) {
				RunUpDiagonal(queryChars, subjectChars, rowMinima, colMinima, minRow - 1, minCol - 1);
			}
		}

		//	Cuts diagonally across the distance matrix in the direction of decreasing 
		//	query index and increasing subject index.
		// This is suited for the case where m < n
		void DoSliceHighLow(
			const Symbol* queryChars,
			const Symbol* subjectChars,
			vector<Distance> & rowMinima,
			vector<Distance> & colMinima,
			vector<bool> & diagColsProcessed,
			const size_t m,
			const size_t n,
			const size_t r,
			const size_t c
		) {
			Distance minDist = MAX_DIST;
			int minRow = 0;
			int minCol = 0;

			for (int row = (int)r, col = (int)c; (size_t)row < m && col >= 0; col--, row++) {
				Distance distance = GetKmerDistance( queryChars + row, subjectChars + col, kmerLength);

				if (distance < minDist) {
					minDist = distance;
					minRow = row;
					minCol = col;
				}

				if (distance < rowMinima[row]) rowMinima[row] = distance;
				if (distance < colMinima[col]) colMinima[col] = distance;

				if (col == 0) continue;

				distance = GetKmerDistance(queryChars + row, subjectChars + col - 1, kmerLength);

				if (distance < minDist) {
					minDist = distance;
					minRow = row;
					minCol = col - 1;
				}

				if (distance < rowMinima[row]) rowMinima[row] = distance;
				if (distance < colMinima[col - 1]) colMinima[col - 1] = distance;
			}

			if (diagColsProcessed[m + minCol - minRow]) return;

			diagColsProcessed[m + minCol - minRow] = true;

			if ((size_t)minRow < m - 1 && (size_t)minCol < n - 1) {
				RunDownDiagonal(queryChars, subjectChars, rowMinima, colMinima, m, n, minRow + 1, minCol + 1);
			}

			if (minRow > 0 && minCol > 0) {
				RunUpDiagonal(queryChars, subjectChars, rowMinima, colMinima, minRow - 1, minCol - 1);
			}
		}

		void RunDownDiagonal(
			const Symbol* queryChars,
			const Symbol* subjectChars,
			vector<Distance> & rowMinima,
			vector<Distance> & colMinima,
			const size_t m,
			const size_t n,
			const size_t r,
			const int c
		) {
			Distance buffer[kmerLength];
			const auto * a = subjectChars + c;
			const auto * b = queryChars + r;
			Distance distance = 0;

			size_t diagLength = std::min(m - r, n - c);

			// Prime the circular buffer with the first kmer in the query
			for (size_t t = 0; t < kmerLength; t++, a++, b++) {
				Distance currentTerm = distanceLookup[(*a).value][(*b).value];
				distance += currentTerm;
				buffer[t] = currentTerm;
			}

			if (distance < rowMinima[r]) rowMinima[r] = distance;
			if (distance < colMinima[c]) colMinima[c] = distance;

			for (size_t offset = 1, buffptr = 0;
				offset < diagLength;
				a++, b++, offset++, buffptr++
				) {
				if (buffptr >= kmerLength) {
					buffptr = 0;
				}

				distance -= buffer[buffptr];
				Distance currentTerm = distanceLookup[( *a ).value][( *b ).value];
				buffer[buffptr] = currentTerm;
				distance += currentTerm;

				if (distance < rowMinima[r + offset]) rowMinima[r + offset] = distance;
				if (distance < colMinima[c + offset]) colMinima[c + offset] = distance;
			}
		}

		void RunUpDiagonal(
			const Symbol* queryChars,
			const Symbol* subjectChars,
			vector<Distance> & rowMinima,
			vector<Distance> & colMinima,
			const int r,
			const int c
		) {
			Distance buffer[kmerLength];
			const Symbol* a = subjectChars + c + kmerLength - 1;
			const Symbol* b = queryChars + r + kmerLength - 1;
			Distance distance = 0;

			// Prime the circular buffer with the first kmer in the query
			for (int t = kmerLength - 1; t >= 0; t--, a--, b--) {
				Distance currentTerm = distanceLookup[a->value][b->value];
				distance += currentTerm;
				buffer[t] = currentTerm;
			}

			if (distance < rowMinima[r]) rowMinima[r] = distance;
			if (distance < colMinima[c]) colMinima[c] = distance;

			int offsetEnd = std::max(-r - 1, -c - 1);

			for (int offset = -1, buffptr = kmerLength - 1;
				offset > offsetEnd;
				a--, b--, offset--, buffptr--
				) {
				if (buffptr < 0) {
					buffptr = kmerLength - 1;
				}

				distance -= buffer[buffptr];
				Distance currentTerm = distanceLookup[(*a).value][(*b).value];
				buffer[buffptr] = currentTerm;
				distance += currentTerm;

				if (distance < rowMinima[r + offset]) rowMinima[r + offset] = distance;
				if (distance < colMinima[c + offset]) colMinima[c + offset] = distance;
			}
		}

#if 0
	protected:
		virtual double GetSequenceDistance(
			Distance * rowMinima,
			Distance * colMinima,
			size_t queryKmerCount,
			size_t subjectKmerCount
		) {
			int rowTotal = 0;
			Distance worst = 0;

			for (uint i = 0; i < queryKmerCount; i++) {
				rowTotal += rowMinima[i];

				if (rowMinima[i] > worst) worst = rowMinima[i];
			}

			int colTotal = 0;

			for (uint i = 0; i < subjectKmerCount; i++) {
				colTotal += colMinima[i];

				if (colMinima[i] > worst) worst = colMinima[i];
			}

			if (fragMode == FragmentAggregationMode::HausdorffAverageAverage()) {
				return ((double)rowTotal / queryKmerCount + (double)colTotal / subjectKmerCount) / 2;
			}
			else if (fragMode == FragmentAggregationMode::HausdorffAverage()) {
				return max((double)rowTotal / queryKmerCount, (double)colTotal / subjectKmerCount);
			}
			else if (fragMode == FragmentAggregationMode::Hausdorff()) {
				return worst;
			}
			else {
				throw NotImplementedException(FileAndLine);
			}
		}

#endif // 0
	};
}
