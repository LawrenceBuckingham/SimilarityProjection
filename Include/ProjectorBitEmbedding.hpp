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
#include "CharMap.hpp"

#include <string>
#include <ctime>
#include <cfloat>

#undef max
#undef min

namespace QutBio {

	class ProjectorBitEmbedding : public Projector {
	protected:
		const CharMap & queryCharMap;
		const CharMap & subjectCharMap;

	public:
		ProjectorBitEmbedding(
			SimilarityMatrix * matrix,
			uint kmerLength,
			const FragmentAggregationMode * fragMode
			) :
			Projector(matrix, kmerLength, fragMode),
			queryCharMap(CharMap::Blosum62QueryEncoding()),
			subjectCharMap(CharMap::Blosum62SubjectEncoding()) {
			// TODO: Change this to general matrix.
		}

		virtual ~ProjectorBitEmbedding() {}

		// overrides the inherited version.
		double ComputeDistance(
			EncodedFastaSequence * querySeq,
			EncodedFastaSequence * subjectSeq
			) {
			auto & queryStr = querySeq->Sequence();
			size_t queryKmerCount = queryStr.size() - kmerLength + 1;

			if ( querySeq->Embedding().size() == 0 ) {
				querySeq->SetEmbedding(queryCharMap);
			}

			auto & subjectStr = subjectSeq->Sequence();
			size_t subjectKmerCount = subjectStr.size() - kmerLength + 1;

			if ( subjectSeq->Embedding().size() == 0 ) {
				subjectSeq->SetEmbedding(subjectCharMap);
			}

			auto rowMinima = (Distance *) alloca(queryKmerCount * sizeof(Distance));
			auto colMinima = (Distance *) alloca(subjectKmerCount * sizeof(Distance));

			auto queryChars = (uint64_t*) alloca(queryStr.size() * sizeof(uint64_t));
			auto subjectChars = (uint64_t*) alloca(subjectStr.size() * sizeof(uint64_t));

			memcpy(queryChars, querySeq->Embedding().data(), queryStr.size() * sizeof(uint64_t));
			memcpy(subjectChars, subjectSeq->Embedding().data(), subjectStr.size() * sizeof(uint64_t));

			ComputeDistanceMatrix(
				queryChars, subjectChars, rowMinima, colMinima, queryKmerCount, subjectKmerCount
				);

			double distance = GetSequenceDistance(rowMinima, colMinima, queryKmerCount, subjectKmerCount);

			return distance;
		}

		/// <summary>
		///	Given strings s and t, having length m and n respectively, and a similarity matrix 
		///	reformatted as a 128 by 128 array of integers, populates a vector contain a kmer mutual
		///	distance table.
		/// </summary>
		///	<param name="queryBytes"></param>

	protected:
		// Override UpdateKmerDistanceCache.
		void ComputeDistanceMatrix(
			const uint64_t * queryChars,
			const uint64_t * subjectChars,
			Distance * rowMinima,
			Distance * colMinima,
			size_t queryKmerCount,
			size_t subjectKmerCount
			) {
			const size_t m = queryKmerCount;
			const size_t n = subjectKmerCount;

			memset(rowMinima, 127, m * sizeof(rowMinima[0]));
			memset(colMinima, 127, n * sizeof(colMinima[0]));

			for ( size_t r = 0; r < m; r++ ) {
				const int c_upper = r == 0 ? (int) n : 1;

				//	Do the top-right part of the rectangle
				for ( int c = 0; c < c_upper; c++ ) {
					Distance buffer[1000];
					const uint64_t * a = subjectChars + c;
					const uint64_t * b = queryChars + r;
					Distance distance = 0;

					size_t diagLength = std::min(m - r, n - c);

					// Prime the circular buffer with the first kmer in the query
					for ( size_t t = 0; t < kmerLength; t++, a++, b++ ) {
						Distance currentTerm = (Distance) POPCOUNT((*a) ^ (*b)) >> 1;
						distance += currentTerm;
						buffer[t] = currentTerm;
					}

					if ( distance < rowMinima[r] ) rowMinima[r] = distance;
					if ( distance < colMinima[c] ) colMinima[c] = distance;

					for ( size_t offset = 1, buffptr = 0;
						offset < diagLength;
						a++, b++, offset++, buffptr++
						) {
						if ( buffptr >= kmerLength ) {
							buffptr = 0;
						}

						distance -= buffer[buffptr];
						Distance currentTerm = (Distance) POPCOUNT((*a) ^ (*b)) >> 1;
						buffer[buffptr] = currentTerm;
						distance += currentTerm;

						if ( distance < rowMinima[r + offset] ) rowMinima[r + offset] = distance;
						if ( distance < colMinima[c + offset] ) colMinima[c + offset] = distance;
					}
				}
			}
		}

	};
}
