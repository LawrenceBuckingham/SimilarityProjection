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
#include "AveragePrecision.hpp"
#include "Util.hpp"

#include "KmerDistanceCalculatorDP.hpp"

namespace QutBio {

	template <typename DistanceFunction, typename KmerType>
	class KmerSequenceRankerDP : public KmerDistanceCalculatorDP<DistanceFunction, KmerType> {
		using BaseClass = KmerDistanceCalculatorDP<DistanceFunction, KmerType>;

	public:
		KmerSequenceRankerDP(
			/**/	SimilarityMatrix * matrix,
			/**/	FragmentAggregationMode * kmerMode,
			/**/	FragmentAggregationMode * fragMode,
			/**/	Alphabet * alphabet,
			/**/	uint fragLength,
			/**/	uint kmerLength,
			/**/	bool pushKmerDistances,
			/**/	DistanceFunction & distanceFunction,
			/**/	int skip,
			/**/	int maxRecords
			) :
			KmerDistanceCalculatorDP<DistanceFunction, KmerType>(
			/**/	matrix,
			/**/	kmerMode,
			/**/	fragMode,
			/**/	alphabet,
			/**/	fragLength,
			/**/	kmerLength,
			/**/	pushKmerDistances,
			/**/	distanceFunction,
			/**/	skip,
			/**/	maxRecords
			) {
		}

		virtual ~KmerSequenceRankerDP() {
		}

#if defined(WANT_DIAGNOSTIC_STREAM)
		void SetDiagnosticStream( ostream * diagnosticStream ) {
			this->diagnosticStream = diagnosticStream;
		}
#endif

		void Setup( void ) {}

		void PreProcess( EncodedFastaSequence & querySeq, int queryIdx ) {
			(*this->rankings)[0].clear();
		}

		void Process( EncodedFastaSequence & querySeq, size_t queryIdx, EncodedFastaSequence & subjectSeq, size_t subjectIdx, double distance ) {
			(*this->rankings)[0].emplace_back(
				querySeq.IdStr(),
				subjectSeq.IdStr(),
				distance,
				0, // rank
				0  // Hits
			);
		}

		void PostProcess( EncodedFastaSequence & querySeq, int queryIdx ) {
		}

		void Cleanup( void ) {
			BaseClass::RunComplete();
		}

		void SetQueryComplete(function<void(vector<Ranking>& rankings)> & value) {
			BaseClass::SetQueryComplete(value);
		}

		void RunJob(
				vector<EncodedFastaSequence *> & query,
				vector<EncodedFastaSequence *> & db
			) {
			BaseClass::RunJob(query, db);
		}
	};
}
