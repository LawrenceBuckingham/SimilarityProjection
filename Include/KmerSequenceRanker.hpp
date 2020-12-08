#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <ctime>
#include <cfloat>
#include <limits>

#include "Alphabet.hpp"
#include "Console.hpp"
#include "Delegates.hpp"
#include "EncodedKmer.hpp"
#include "FastaSequence.hpp"
#include "Fragment.hpp"
#include "KmerSequenceRanker_Params.hpp"
#include "Ranking.hpp"
#include "AveragePrecision.hpp"
#include "Util.hpp"

#include "KmerDistanceCalculator.hpp"

namespace QutBio {

	template <typename DistanceFunction, typename KmerType>
	class KmerSequenceRanker : public KmerDistanceCalculator<DistanceFunction, KmerType> {

		using Codebook = KmerCodebook<DistanceFunction, KmerType>;

	public:
		KmerSequenceRanker(
			/**/	SimilarityMatrix * matrix,
			/**/	FragmentAggregationMode * kmerMode,
			/**/	FragmentAggregationMode * fragMode,
			/**/	Alphabet * alphabet,
			/**/	uint fragLength,
			/**/	uint kmerLength,
			/**/	Distance thresholdDistance,
			/**/	Distance defaultDistance,
			/**/	bool pushKmerDistances,
			/**/	DistanceFunction & distanceFunction,
			/**/	int skip,
			/**/	int maxRecords
			) :
			KmerDistanceCalculator<DistanceFunction, KmerType>(
			/**/	matrix,
			/**/	kmerMode,
			/**/	fragMode,
			/**/	alphabet,
			/**/	fragLength,
			/**/	kmerLength,
			/**/	thresholdDistance,
			/**/	defaultDistance,
			/**/	pushKmerDistances,
			/**/	distanceFunction,
			/**/	skip,
			/**/	maxRecords
			) {
		}

		virtual ~KmerSequenceRanker() {
		}

		void SetCodebook( Codebook * codebook ) {
			this->codebook = codebook;
		}

#if defined(WANT_DIAGNOSTIC_STREAM)
		void SetDiagnosticStream( ostream * diagnosticStream ) {
			this->diagnosticStream = diagnosticStream;
		}
#endif

		void setup( void ) {}

		void preProcess( EncodedFastaSequence & querySeq, int queryIdx ) {
			(*KmerDistanceCalculator<DistanceFunction, KmerType>::rankings)[0].clear();
		}

		void process( EncodedFastaSequence & querySeq, size_t queryIdx, EncodedFastaSequence & subjectSeq, size_t subjectIdx, double distance ) {
			(*KmerDistanceCalculator<DistanceFunction, KmerType>::rankings)[0].emplace_back(
				querySeq.IdStr(),
				subjectSeq.IdStr(),
				distance,
				0, // rank
				0  // Hits
			);
		}

		void postProcess( EncodedFastaSequence & querySeq, int queryIdx ) {
		}

		void cleanup( void ) {
			KmerDistanceCalculator<DistanceFunction, KmerType>::runComplete();
		}

		void SetQueryComplete(function<void(vector<Ranking>& rankings)> & value) {
			KmerDistanceCalculator<DistanceFunction, KmerType>::SetQueryComplete(value);
		}

		void RunJob(
				vector<EncodedFastaSequence *> & query,
				vector<EncodedFastaSequence *> & db
			) {
			KmerDistanceCalculator<DistanceFunction, KmerType>::RunJob(query, db);
		}
	};
}
