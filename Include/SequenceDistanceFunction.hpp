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
#include <KmerClusterPrototype.hpp>

#include <string>
#include <ctime>
#include <cfloat>

namespace QutBio {

	class SequenceDistanceFunction {
	protected:
		Distance distanceLookup[128][128];
		uint kmerLength;

	public:

		SequenceDistanceFunction(
			SimilarityMatrix * matrix,
			uint kmerLength
		) :
			kmerLength(kmerLength)
			//
		{
			matrix->PopulateDistanceTable(distanceLookup);
		}

		virtual ~SequenceDistanceFunction() {}

		virtual double ComputeDistance(
			FastaSequence & querySeq,
			FastaSequence & subjectSeq
		) = 0;

		Distance GetKmerDistance(
			const Symbol* queryBytes,
			const Symbol* subjectBytes,
			size_t kmerLength
		) {
			Distance distance = 0;

			for (size_t t = 0; t < kmerLength; t++) {
				distance += distanceLookup[queryBytes[t].value][subjectBytes[t].value];
			}

			return distance;
		}
	};
}
