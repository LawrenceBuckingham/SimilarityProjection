#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "SimilarityMatrix.hpp"

namespace QutBio {
	class DiagonalGenerator  {
	public:
		/**
		 *	<summary>
		 *		Trying to ease the way towards a fully flexible but manageable framework 
		 *		that can deal with bulk nearby neighbours from the entire database.
		 *		The idea is to feed each distance that we discover forward to the consumer
		 *		rather than pre-computing all distances and then looking back.
		 *	<para>
		 *		Not all arguments will be useful for all rankers, and at the moment I am mainly 
		 *		interested in the Projector and BestOfBest+Fragmentation rankers.
		 *	</para>
		 *	</summary>
		 */
		typedef void( *DistanceProcessor )( 
			void * processingObject, 
			size_t queryPos, 
			size_t subjectPos, 
			Distance dist 
		);

		void GenerateDistances(
			// EncodedFastaSequence * querySeq,
			// EncodedFastaSequence * subjectSeq,
			const Symbol* queryChars_,
			const Symbol* subjectChars_,
			size_t kmerLength,
			size_t queryKmerCount,
			size_t subjectKmerCount,
			const Distance distanceLookup[128][128],
			void * processingObject,
			DistanceProcessor process
			) {
			const size_t m = queryKmerCount;
			const size_t n = subjectKmerCount;
			const uint8_t * queryChars = (uint8_t *)queryChars_;
			const uint8_t * subjectChars = (uint8_t *)subjectChars_;

			for ( size_t r = 0; r < m; r++ ) {
				const int c_upper = r == 0 ? (int) n : 1;

				//	Do the top-right part of the rectangle
				for ( int c = 0; c < c_upper; c++ ) {
					Distance buffer[1000];
					const uint8_t * a = subjectChars + c;
					const uint8_t * b = queryChars + r;
					Distance distance = 0;

					size_t diagLength = std::min( m - r, n - c );

					// Prime the circular buffer with the first kmer in the query
					for ( size_t t = 0; t < kmerLength; t++, a++, b++ ) {
						Distance currentTerm = distanceLookup[*a][*b];
						distance += currentTerm;
						buffer[t] = currentTerm;
					}

					process( processingObject, r, c, distance );

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

						process( processingObject, r + offset, c + offset, distance );
					}
				}
			}
		}


	};
}
