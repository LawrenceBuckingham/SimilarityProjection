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

namespace QutBio {

	class BestHitDetector: SequenceDistanceFunction {
	protected:
		/// <summary>
		///		Row and column minima
		///	</summary>
		Distance bestDistance;

		// Query sequence data.
		EncodedFastaSequence * querySeq;
		const byte * queryChars;
		size_t queryIdx;
		size_t queryKmerCount;

		// Subject (DB) sequence data.
		EncodedFastaSequence * subjectSeq;
		const byte * subjectChars;
		size_t subjectIdx;
		size_t subjectKmerCount;

	public:
		BestHitDetector(
			SimilarityMatrix * matrix,
			uint kmerLength
			) : SequenceDistanceFunction( matrix, kmerLength ) {}

		virtual ~BestHitDetector() {}

		double ComputeDistance(
			EncodedFastaSequence * querySeq,
			EncodedFastaSequence * subjectSeq
			) {
			queryChars = querySeq->Sequence().data();
			queryKmerCount = querySeq->Sequence().size() - kmerLength + 1;

			subjectChars = subjectSeq->Sequence().data();
			subjectKmerCount = subjectSeq->Sequence().size() - kmerLength + 1;

			ComputeDistanceMatrix();

			return bestDistance;
		}

	private:

		/// <summary>
		///	Given strings s and t, having length m and n respectively, and a similarity matrix 
		///	reformatted as a 128 by 128 array of integers, populates a vector contain a kmer mutual
		///	distance table.
		/// </summary>
		///	<param name="queryBytes"></param>

		void ComputeDistanceMatrix() {
			const size_t m = queryKmerCount;
			const size_t n = subjectKmerCount;

			bestDistance = numeric_limits<Distance>::max();

			for ( size_t r = 0; r < m; r++ ) {
				const int c_upper = r == 0 ? (int) n : 1;

				//	Do the top-right part of the rectangle
				for ( int c = 0; c < c_upper; c++ ) {
					Distance buffer[1000];
					const uint8_t * a = (uint8_t *) subjectChars + c;
					const uint8_t * b = (uint8_t *) queryChars + r;
					Distance distance = 0;

					size_t diagLength = m - r;

					if ( n - c < diagLength ) {
						diagLength = n - c;
					}

					// Prime the circular buffer with the first kmer in the query
					for ( size_t t = 0; t < kmerLength; t++, a++, b++ ) {
						Distance currentTerm = distanceLookup[*a][*b];
						distance += currentTerm;
						buffer[t] = currentTerm;
					}

					if ( distance < bestDistance ) bestDistance = distance;

					for ( size_t offset = 1, buffptr = 0;
						offset < diagLength;
						a++, b++, offset++, buffptr++
						) {
						if ( buffptr >= kmerLength ) {
							buffptr = 0;
						}

						distance -= buffer[buffptr];
						Distance currentTerm = distanceLookup[*a][*b];
						buffer[buffptr] = currentTerm;
						distance += currentTerm;

						if ( distance < bestDistance ) bestDistance = distance;
					}
				}
			}
		}

	};
}
