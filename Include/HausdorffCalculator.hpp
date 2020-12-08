#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "Alphabet.hpp"
#include "SimilarityMatrix.hpp"
#include "FastaSequence.hpp"
#include "FragDistProcessor.hpp"
#include "Fragment.hpp"
#include "FreeList.hpp"
#include "SequenceDistanceFunction.hpp"
#include "Util.hpp"
#include "DiagonalGenerator.hpp"
#include "FragmentAggregationMode.hpp"

#include <string>
#include <ctime>
#include <cfloat>

#undef max
#undef min

namespace QutBio {

	class HausdorffCalculator : SequenceDistanceFunction {
	public:

		typedef Distance(*FragmentDistance)(
			HausdorffCalculator & This,

			size_t qStart,
			size_t qEnd,

			size_t sStart,
			size_t sEnd
			);

		typedef double(*CollectionDistance) (
			HausdorffCalculator & This,

			size_t queryKmerCount,
			size_t queryFragCount,

			size_t subjectKmerCount,
			size_t subjectFragCount
			);

#if WANT_DIAGNOSTIC_STREAM
	public:
#else
	protected:
#endif
		/// <summary>
		///		kmer distance cache with max(queryKmerCount) rows and max(subjectKmerCount) columns.
		///	</summary>
		vector<Distance> kmerDistCache_;
		MatrixView<Distance> kmerDistCache;
		vector<Distance> rowMinima;
		vector<Distance> colMinima;
		FragmentAggregationMode * kmerMode;
		FragmentAggregationMode * fragMode;
		uint kmerLength;
		uint   fragmentLength;
		size_t maxQueryLength;
		size_t maxSubjectLength;
		const Symbol* queryBytes;
		const Symbol* subjectBytes;

#if defined(USE_BIG_CACHE)
		int vocabSize;
		size_t charsPerWord;
		KmerDistanceCache & cachedCalculator;
		KmerDistanceFunction kmerCodeDistance;
#endif // defined(USE_BIG_CACHE)
		FragmentDistance fragmentDistance;
		CollectionDistance collectionDistance;

		Distance thresholdDistance = BAD_DIST;
		Distance defaultDistance = MAX_DIST;

	public:
		HausdorffCalculator(
			SimilarityMatrix * matrix,
			uint kmerLength,
			FragmentAggregationMode * kmerMode,
			FragmentAggregationMode * fragMode,
			Alphabet * alphabet,
			uint fragLength_,
			size_t maxQueryLength_,
			size_t maxSubjectLength_
		) :
			SequenceDistanceFunction(matrix, kmerLength),
			kmerDistCache_(maxQueryLength_ * maxSubjectLength_),
			kmerDistCache(kmerDistCache_.data(), maxQueryLength_, maxSubjectLength_),
			rowMinima(Fragment::GetCount(maxQueryLength_, fragLength_)),
			colMinima(Fragment::GetCount(maxSubjectLength_, fragLength_)),
			kmerMode(kmerMode),
			fragMode(fragMode),
			kmerLength(kmerLength),
			fragmentLength(fragLength_),
			maxQueryLength(maxQueryLength_),
			maxSubjectLength(maxSubjectLength_)
			//
		{

			FragmentDistance fragmentDistances[] = {
				KmerBestOfBest,
				KmerHausdorff,
				KmerHausdorffAverage,
				KmerHausdorffAverageAverage,
				KmerSlice,
				KmerSliceVertical,
				KmerSliceNoFollow,
				KmerSliceVerticalNoFollow,
			};

			fragmentDistance = fragmentDistances[kmerMode->Value()];

			CollectionDistance collectionDistances[] = {
				FragBestOfBest,
				FragHausdorff,
				FragHausdorffAverage,
				FragHausdorffAverageAverage
			};

			collectionDistance = collectionDistances[fragMode->Value()];
		}

		virtual ~HausdorffCalculator() {}

		void SetThreshold(Distance thresholdDistance, Distance defaultDistance) {
			assert_false(IS_BAD_SIM(thresholdDistance));
			assert_false(IS_BAD_SIM(defaultDistance));
			assert_true(thresholdDistance <= defaultDistance);

			this->thresholdDistance = thresholdDistance;
			this->defaultDistance = defaultDistance;
		}

		double ThresholdDistance() { return thresholdDistance; }

		double DefaultDistance() { return defaultDistance; }

		size_t    queryFragCount;
		size_t subjectFragCount;

		size_t    queryKmerCount;
		size_t subjectKmerCount;

		/**
		 *	Populates the kmer distance matrix, row-minima, and column-minima,
		 *	using the fragment distance function to combine kmer distances from
		 *	the matrix to fragment distances, and using the collection distance
		 *	to combine fragment distances, yielding row- and column-minima and
		 *	the final combined sequence distance.
		 */
		double ComputeDistance(
			FastaSequence & querySeq,
			FastaSequence & subjectSeq
		) {
			auto &queryStr = querySeq.Sequence();
			queryKmerCount = querySeq.KmerCount(kmerLength);

			auto & subjectStr = subjectSeq.Sequence();
			subjectKmerCount = subjectSeq.KmerCount(kmerLength);

			queryFragCount = Fragment::GetCount(queryKmerCount, fragmentLength);
			subjectFragCount = Fragment::GetCount(subjectKmerCount, fragmentLength);

			rowMinima.resize(queryFragCount);
			colMinima.resize(subjectFragCount);

			std::fill(rowMinima.begin(), rowMinima.end(), defaultDistance);
			std::fill(colMinima.begin(), colMinima.end(), defaultDistance);

			kmerDistCache.Reinterpret(queryKmerCount, subjectKmerCount);

			ComputeDistanceMatrix(queryStr.data(), subjectStr.data(), queryKmerCount, subjectKmerCount, kmerLength, distanceLookup, kmerDistCache);

			double distance = collectionDistance(
				*this,
				queryKmerCount, queryFragCount,
				subjectKmerCount, subjectFragCount
			);

			return distance;
		}

		/// <summary>
		/// Gets a reference to the pairwise kmer distance matrix
		/// </summary>
		MatrixView<Distance> & GetKmerDistances() {
			return kmerDistCache;
		}

		/// <summary>
		/// Gets a reference to the row minima vector
		/// </summary>
		vector<Distance> & GetRowMinima() {
			return rowMinima;
		}

		/// <summary>
		/// Gets a reference to the column minima vector
		/// </summary>
		vector<Distance> & GetColMinima() {
			return colMinima;
		}

		/// <summary>
		///	Given strings s and t, having length m and n respectively, and a similarity matrix 
		///	reformatted as a 128 by 128 array of integers, populates a vector contain a kmer mutual
		///	distance table.
		/// </summary>
		///	<param name="queryBytes"></param>

		static void ComputeDistanceMatrix(
			const Symbol* queryChars_,
			const Symbol* subjectChars_,
			size_t queryKmerCount,
			size_t subjectKmerCount,
			size_t kmerLength,
			Distance lookup[128][128],
			MatrixView<Distance> &kmerDistCache
		) {
			const size_t m = queryKmerCount;
			const size_t n = subjectKmerCount;

			if (m == 0 || n == 0) {
				cerr << "Sequence with 0 kmers encountered.\n"
					"query: " << queryChars_ << "\n"
					"reference: " << subjectChars_ << "\n";
				return;
			}

			const uint8_t * queryChars = (uint8_t *)queryChars_;
			const uint8_t * subjectChars = (uint8_t *)subjectChars_;

			for (int r = 0; r < (int)m; r++) {
				const int c_upper = r == 0 ? (int)n : 1;

				//	Do the top-right part of the rectangle
				for (int c = 0; c < c_upper; c++) {
					Distance buffer[1000];
					const uint8_t * a = subjectChars + c;
					const uint8_t * b = queryChars + r;
					Distance distance = 0;

					size_t diagLength = std::min(m - r, n - c);

					// Prime the circular buffer with the first kmer in the query
					for (size_t t = 0; t < kmerLength; t++, a++, b++) {
						Distance currentTerm = lookup[*a][*b];
						distance += currentTerm;
						buffer[t] = currentTerm;
					}

					kmerDistCache(r, c) = distance;

					for (size_t offset = 1, buffptr = 0; offset < diagLength; a++, b++, offset++, buffptr++) {
						if (buffptr >= kmerLength) {
							buffptr = 0;
						}

						distance -= buffer[buffptr];
						Distance currentTerm = lookup[*a][*b];
						buffer[buffptr] = currentTerm;
						distance += currentTerm;

						kmerDistCache(r + offset, c + offset) = distance;
					}
				}
			}
		}

	public:

		///	<summary>
		///		Calculates the threshold Hausdorff average (bidirectional) 
		///		distance between two fragments using the maximum of the two
		///		one way distances.
		///	</summary>
		static double FragHausdorffAverage(
			HausdorffCalculator & This,

			size_t queryKmerCount,
			size_t queryFragCount,

			size_t subjectKmerCount,
			size_t subjectFragCount
		) {
			double totXY = 0;
			int obsXY = 0;

			double totYX = 0;
			int obsYX = 0;

			// Computes distance between fragment at qFragIdx in query and fragment at sFragIdx in subject.
			// Indices run from 0 to queryFragCount in the query, and from 0 to sFragCount in the subject.
			auto processCell = [&](uint qFragIdx, uint qStart, uint qEnd, uint sFragIdx, uint sStart, uint sEnd) {
				auto distance = This.fragmentDistance(This, qStart, qEnd, sStart, sEnd);
				if (distance < This.rowMinima[qFragIdx]) This.rowMinima[qFragIdx] = distance;
				if (distance < This.colMinima[sFragIdx]) This.colMinima[sFragIdx] = distance;
			};

			// Updates the running total after each fragment in the query has been compared with all
			// fragments in the subject.
			auto processEndOfRow = [&](uint qFragIdx) {
				obsXY++;
				totXY += This.rowMinima[qFragIdx];
			};

			// Does the map-reduce operation.
			Fragment::PartitionSequencePair(This.fragmentLength,
				queryKmerCount, queryFragCount,
				subjectKmerCount, subjectFragCount,
				processCell, processEndOfRow
			);

			// Compute one-way distance from y to x, using cache distances.
			for (uint j = 0; j < subjectFragCount; j++) {
				obsYX++;
				totYX += This.colMinima[j];
			}

			double avgXY = totXY / obsXY;
			double avgYX = totYX / obsYX;
			double result = avgXY > avgYX ? avgXY : avgYX;
			return result;
		}

		///	<summary>
		///		Calculates the threshold Hausdorff average (bidirectional) 
		///		distance between two fragments using the average of the two 
		///		one-way fragment distances.
		///	</summary>
		static double FragHausdorffAverageAverage(
			HausdorffCalculator & This,

			size_t queryKmerCount,
			size_t queryFragCount,

			size_t subjectKmerCount,
			size_t subjectFragCount
		) {
			double totXY = 0;
			int obsXY = 0;

			double totYX = 0;
			int obsYX = 0;


			// Computes distance between fragment at qFragIdx in query and fragment at sFragIdx in subject.
			// Indices run from 0 to queryFragCount in the query, and from 0 to sFragCount in the subject.
			auto processCell = [&](uint qFragIdx, uint qStart, uint qEnd, uint sFragIdx, uint sStart, uint sEnd) {
				auto distance = This.fragmentDistance(This, qStart, qEnd, sStart, sEnd);
				if (distance < This.rowMinima[qFragIdx]) This.rowMinima[qFragIdx] = distance;
				if (distance < This.colMinima[sFragIdx]) This.colMinima[sFragIdx] = distance;
			};

			// Updates the running total after each fragment in the query has been compared with all
			// fragments in the subject.
			auto processEndOfRow = [&](uint qFragIdx) {
				obsXY++;
				totXY += This.rowMinima[qFragIdx];
			};

			// Does the map-reduce operation.
			Fragment::PartitionSequencePair(This.fragmentLength,
				queryKmerCount, queryFragCount,
				subjectKmerCount, subjectFragCount,
				processCell, processEndOfRow
			);

			// Compute one-way distance from y to x, using cache distances.
			for (uint j = 0; j < subjectFragCount; j++) {
				obsYX++;
				totYX += This.colMinima[j];
			}

			double avgXY = totXY / obsXY;
			double avgYX = totYX / obsYX;

			double result = (avgXY + avgYX) / 2;
			return result;
		}

		static double FragHausdorff(
			HausdorffCalculator & This,

			size_t queryKmerCount,
			size_t queryFragCount,

			size_t subjectKmerCount,
			size_t subjectFragCount
		) {
			double maxXY = -DBL_MAX;
			double maxYX = -DBL_MAX;

			double qStepSize = Fragment::GetRealStepSize(queryKmerCount, This.fragmentLength, queryFragCount);
			double sStepSize = Fragment::GetRealStepSize(subjectKmerCount, This.fragmentLength, subjectFragCount);

			// Compute one-way distance from x to y, and also cache distances.

			uint qStart = 0;

			for (int i = 0; i < (int)queryFragCount; i++) {
				uint qEnd = Fragment::GetFragmentStart(i + 1, qStepSize, queryKmerCount);

				uint sStart = 0;

				for (int j = 0; j < (int)subjectFragCount; j++) {
					uint sEnd = Fragment::GetFragmentStart(j + 1, sStepSize, subjectKmerCount);

					auto distance = This.fragmentDistance(This, qStart, qEnd, sStart, sEnd);

					if (distance < This.rowMinima[i]) This.rowMinima[i] = distance;
					if (distance < This.colMinima[j]) This.colMinima[j] = distance;

					sStart = sEnd;
				}

				qStart = qEnd;
			}

			// Compute one-way distance from x to y, using cache distances.
			for (uint i = 0; i < queryFragCount; i++) {
				if (This.rowMinima[i] > maxXY) {
					maxXY = This.rowMinima[i];
				}
			}

			// Compute one-way distance from y to x, using cache distances.
			for (uint j = 0; j < subjectFragCount; j++) {
				if (This.colMinima[j] > maxYX) {
					maxYX = This.colMinima[j];
				}
			}

			return maxXY > maxYX ? maxXY : maxYX;
		}

		static double FragBestOfBest(
			HausdorffCalculator & This,

			size_t queryKmerCount,
			size_t queryFragCount,

			size_t subjectKmerCount,
			size_t subjectFragCount
		) {
			double minDist = numeric_limits<double>::max();

			double qStepSize = Fragment::GetRealStepSize(queryKmerCount, This.fragmentLength, queryFragCount);
			double sStepSize = Fragment::GetRealStepSize(subjectKmerCount, This.fragmentLength, subjectFragCount);

			// Compute one-way distance from x to y, and also cache distances.

			uint qStart = 0;

			for (uint i = 0; i < queryFragCount; i++) {
				uint qEnd = Fragment::GetFragmentStart(i + 1, qStepSize, queryKmerCount);

				uint sStart = 0;

				for (uint j = 0; j < subjectFragCount; j++) {
					uint sEnd = Fragment::GetFragmentStart(j + 1, sStepSize, subjectKmerCount);

					auto distance = This.fragmentDistance(This, qStart, qEnd, sStart, sEnd);

					if (distance < This.rowMinima[i]) This.rowMinima[i] = distance;
					if (distance < This.colMinima[j]) This.colMinima[j] = distance;

					sStart = sEnd;
				}

				qStart = qEnd;
			}

			for (uint i = 0; i < queryFragCount; i++) {
				if (This.rowMinima[i] < minDist) {
					minDist = This.rowMinima[i];
				}
			}

			return minDist;
		}

		/// <summary> Returns the min-of-min kmer distance between two sets of encoded kmers,
		/// </summary>
		/// <param name="qStart">The index of the first kmer in the query fragment.</param>
		/// <param name="qEnd">The index of the first kmer AFTER the end of the query fragment.</param>
		/// <param name="sStart">The index of the first kmer in the reference fragment.</param>
		/// <param name="sEnd">The index of the first kmer AFTER the end of the reference fragment.</param>
		/// <param name="len">The number of kmers in the fragment.</param>
		///	<param name="threshold">Early stopping criterion for summative, non-negative, distances. 
		///		If cumulative distance equals or exceeds threshold, stop calculating and return cumulative 
		///		value at the stopping time.
		///	</param>
		/// <returns></returns>

		static Distance KmerBestOfBest(
			HausdorffCalculator & This,

			size_t qStart,
			size_t qEnd,

			size_t sStart,
			size_t sEnd
		) {
			Distance minDist = numeric_limits<Distance>::max();

			for (size_t i = qStart; i < qEnd; i++) {
				auto row = This.kmerDistCache.Row(i);

				for (size_t j = sStart; j < sEnd; j++) {
					Distance distance = row[j];

					if (distance < minDist) {
						minDist = distance;
					}
				}
			}

			if (!IS_BAD_DIST(This.thresholdDistance)) {
				return minDist <= This.thresholdDistance ? minDist : This.defaultDistance;
			}
			else {
				return minDist;
			}
		}

		static Distance KmerSlice(
			HausdorffCalculator & This,

			size_t qStart,
			size_t qEnd,

			size_t sStart,
			size_t sEnd
		) {
			Distance minDist = numeric_limits<Distance>::max();
			size_t minI = 0;
			size_t minJ = 0;

			for (size_t i = qStart; i < qEnd; i++) {
				size_t j0 = (sEnd - 1) - (i - qStart) * (sEnd - sStart) / (qEnd - qStart);

				for (int j = (int)j0; j >= (int)sStart && j >= (int)j0 - 1; j--) {
					// assert_true((int) sStart <= j && j < (int) sEnd);
					auto distance = This.GetKmerDistance(This.queryBytes+i, This.subjectBytes+j, This.kmerLength);

					if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
						distance = This.defaultDistance;
					}

					if (distance < minDist) {
						minDist = distance;
						minI = i;
						minJ = j;
					}
				}
			}

			for (size_t i = minI + 1, j = minJ + 1; i < qEnd && j < sEnd; i++, j++) {
				Distance distance = This.GetKmerDistance(This.queryBytes + i, This.subjectBytes + j, This.kmerLength);

				if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
					distance = This.defaultDistance;
				}

				if (distance < minDist) {
					minDist = distance;
				}
			}

			for (int i = (int)minI - 1, j = (int)minJ - 1; i >= (int)qStart && j >= (int)sStart; i--, j--) {
				Distance distance = This.GetKmerDistance(This.queryBytes + i, This.subjectBytes + j, This.kmerLength);

				if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
					distance = This.defaultDistance;
				}

				if (distance < minDist) {
					minDist = distance;
				}
			}

			return minDist;
		}

		static Distance KmerSliceVertical(
			HausdorffCalculator & This,

			size_t qStart,
			size_t qEnd,

			size_t sStart,
			size_t sEnd
		) {
			Distance minDist = numeric_limits<Distance>::max();;
			size_t minI = 0;
			size_t minJ = 0;

			for (size_t i = qStart; i < qEnd; i++) {
				size_t j = sStart + (sEnd - sStart) / 2;

				Distance distance = This.GetKmerDistance(This.queryBytes + i, This.subjectBytes + j, This.kmerLength);

				if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
					distance = This.defaultDistance;
				}

				if (distance < minDist) {
					minDist = distance;
					minI = i;
					minJ = j;
				}
			}

			for (size_t i = minI + 1, j = minJ + 1; i < qEnd && j < sEnd; i++, j++) {
				Distance distance = This.GetKmerDistance(This.queryBytes + i, This.subjectBytes + j, This.kmerLength);

				if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
					distance = This.defaultDistance;
				}

				if (distance < minDist) {
					minDist = distance;
				}
			}

			for (int i = (int)minI - 1, j = (int)minJ - 1; i >= (int)qStart && j >= (int)sStart; i--, j--) {
				Distance distance = This.GetKmerDistance(This.queryBytes + i, This.subjectBytes + j, This.kmerLength);

				if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
					distance = This.defaultDistance;
				}

				if (distance < minDist) {
					minDist = distance;
				}
			}

			return minDist;
		}

		static Distance KmerSliceNoFollow(
			HausdorffCalculator & This,

			size_t queryKmerCount,
			size_t qStart,

			size_t subjectKmerCount,
			size_t sStart
		) {
			Distance minDist = numeric_limits<Distance>::max();
			size_t qEnd = min(qStart + This.fragmentLength, queryKmerCount);
			size_t sEnd = min(sStart + This.fragmentLength, subjectKmerCount);

			for (size_t i = qStart; i < qEnd; i++) {
				size_t j0 = (sEnd - 1) - (i - qStart) * (sEnd - sStart) / (qEnd - qStart);

				for (int j = (int)j0; j >= (int)sStart && j >= (int)j0 - 1; j--) {
					assert_true((int)sStart <= j && j < (int)sEnd);

					Distance distance = This.GetKmerDistance(This.queryBytes + i, This.subjectBytes + j, This.kmerLength);

					if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
						distance = This.defaultDistance;
					}

					if (distance < minDist) {
						minDist = distance;
					}
				}
			}

			return minDist;
		}

		static Distance KmerSliceVerticalNoFollow(
			HausdorffCalculator & This,

			size_t queryKmerCount,
			size_t qStart,

			size_t subjectKmerCount,
			size_t sStart
		) {
			Distance minDist = numeric_limits<Distance>::max();
			size_t qEnd = min(qStart + This.fragmentLength, queryKmerCount);
			size_t sEnd = min(sStart + This.fragmentLength, subjectKmerCount);

			for (size_t i = qStart; i < qEnd; i++) {
				size_t j = sStart + (sEnd - sStart) / 2;

				Distance distance = This.GetKmerDistance(This.queryBytes + i, This.subjectBytes + j, This.kmerLength);

				if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
					distance = This.defaultDistance;
				}

				if (distance < minDist) {
					minDist = distance;
				}
			}

			return minDist;
		}

		/// <summary> Returns the max(max-of-min) kmer distance between two sets of encoded kmers,
		/// </summary>
		/// <param name="qStart">The index of the first kmer in the query fragment.</param>
		///	<param name="qEnd">The index of the first kmer AFTER the end of the query fragment.</param>
		/// <param name="sStart">The index of the first kmer in the reference fragment.</param>
		///	<param name="sEnd">The index of the first kmer AFTER the end of the reference fragment.</param>
		/// <param name="len">The number of kmers in the fragment.</param>
		///	<param name="threshold">Early stopping criterion for summative, non-negative, distances. 
		///		If cumulative distance equals or exceeds threshold, stop calculating and return cumulative 
		///		value at the stopping time.
		///	</param>
		/// <returns></returns>

		static Distance KmerHausdorff(
			HausdorffCalculator & This,

			size_t qStart,
			size_t qEnd,

			size_t sStart,
			size_t sEnd
		) {
			Distance maxXY = numeric_limits<short>::min();
			Distance maxYX = numeric_limits<short>::min();

			// Compute one-way distance from x to y, and also cache distances.
			for (size_t i = qStart; i < qEnd; i++) {
				Distance min = numeric_limits<Distance>::max();
				auto row = This.kmerDistCache.Row(i);

				for (size_t j = sStart; j < sEnd; j++) {
					Distance distance = row[j];

					if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
						distance = This.defaultDistance;
					}

					if (distance < min) {
						min = distance;
					}
				}

				if (min > maxXY) {
					maxXY = min;
				}
			}

			// Compute one-way distance from y to x, using cache distances.
			for (size_t j = sStart; j < sEnd; j++) {
				Distance min = numeric_limits<Distance>::max();

				for (size_t i = qStart; i < qEnd; i++) {
					Distance distance = This.kmerDistCache(i, j);

					if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
						distance = This.defaultDistance;
					}

					if (distance < min) {
						min = distance;
					}
				}

				if (min > maxYX) {
					maxYX = min;
				}
			}

			return maxXY > maxYX ? maxXY : maxYX;
		}

		/// <summary> Returns the max(average-of-min) kmer distance between two sets of encoded kmers,
		/// </summary>
		/// <param name="qStart">The index of the first kmer in the query fragment.</param>
		/// <param name="qEnd">The index of the first kmer AFTER the end of the query fragment.</param>
		/// <param name="sStart">The index of the first kmer in the reference fragment.</param>
		/// <param name="sEnd">The index of the first kmer AFTER the end of the reference fragment.</param>
		/// <param name="len">The number of kmers in the fragment.</param>
		/// <returns></returns>

		static Distance KmerHausdorffAverage(
			HausdorffCalculator & This,

			size_t qStart,
			size_t qEnd,

			size_t sStart,
			size_t sEnd
		) {
			int totXY = 0;
			int totYX = 0;

			// Compute one-way distance from x to y, using pre-cached cache distances.
			for (size_t i = qStart; i < qEnd; i++) {
				Distance min = numeric_limits<Distance>::max();
				auto row = This.kmerDistCache.Row(i);

				for (size_t j = sStart; j < sEnd; j++) {
					Distance distance = row[j];

					if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
						distance = This.defaultDistance;
					}

					if (distance < min) {
						min = distance;
					}
				}

				totXY += min;
			}

			// Compute one-way distance from y to x, using cache distances.
			for (size_t j = sStart; j < sEnd; j++) {
				Distance min = numeric_limits<Distance>::max();

				for (size_t i = qStart; i < qEnd; i++) {
					Distance distance = This.kmerDistCache(i, j);

					if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
						distance = This.defaultDistance;
					}

					if (distance < min) {
						min = distance;
					}
				}

				totYX += min;
			}

			double avgXY = (double)totXY / (qEnd - qStart);
			double avgYX = (double)totYX / (sEnd - sStart);

			return avgXY > avgYX ? (Distance)round(avgXY) : (Distance)round(avgYX);
		}

		/// <summary> Returns the average(average-of-min) kmer distance between two sets of encoded kmers,
		/// </summary>
		/// <param name="qStart">The index of the first kmer in the query fragment.</param>
		/// <param name="qEnd">The index of the first kmer AFTER the end of the query fragment.</param>
		/// <param name="sStart">The index of the first kmer in the reference fragment.</param>
		/// <param name="sEnd">The index of the first kmer AFTER the end of the reference fragment.</param>
		/// <returns></returns>

		static Distance KmerHausdorffAverageAverage(
			HausdorffCalculator & This,

			size_t qStart,
			size_t qEnd,

			size_t sStart,
			size_t sEnd
		) {
			int totXY = 0;
			int totYX = 0;

			// Compute one-way distance from x to y, using pre-cached cache distances.
			for (size_t i = qStart; i < qEnd; i++) {
				Distance min = numeric_limits<Distance>::max();
				auto  row = This.kmerDistCache.Row(i);

				for (size_t j = sStart; j < sEnd; j++) {
					Distance distance = row[j];

					if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
						distance = This.defaultDistance;
					}

					if (distance < min) {
						min = distance;
					}
				}

				totXY += min;
			}

			// Compute one-way distance from y to x, using cache distances.
			for (size_t j = sStart; j < sEnd; j++) {
				Distance min = numeric_limits<Distance>::max();

				for (size_t i = qStart; i < qEnd; i++) {
					Distance distance = This.kmerDistCache(i, j);

					if (!IS_BAD_DIST(This.thresholdDistance) && distance > This.thresholdDistance) {
						distance = This.defaultDistance;
					}

					if (distance < min) {
						min = distance;
					}
				}

				totYX += min;
			}

			double avgXY = (double)totXY / (qEnd - qStart);
			double avgYX = (double)totYX / (sEnd - sStart);

			return (Distance)round((avgXY + avgYX) / 2);
		}

	};
}
