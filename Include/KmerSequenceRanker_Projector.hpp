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

#if defined(USE_AMP)
#include "KmerDistanceCalculator_AMP.hpp"
#elif defined(PROJECTOR)
#include "KmerDistanceCalculator_Projector.hpp"
#else
#include "KmerDistanceCalculator.hpp"
#endif

namespace QutBio {

	class KmerSequenceRanker : public KmerDistanceCalculator {

		vector<Ranking> rankings;
		double threshold = HUGE_VAL;
		function<void(vector<Ranking>& rankings)> queryComplete;
		function<void(void)> runComplete;

	public:
		KmerSequenceRanker(
			/**/	DistanceType * dist,
			/**/	uint matrixId,
			/**/	FragmentAggregationMode * kmerMode,
			/**/	FragmentAggregationMode * fragMode,
			/**/	Alphabet * alphabet,
			/**/	uint fragLength,
			/**/	uint interval,
			/**/	uint kmerLength
#if defined(USE_BIG_CACHE)
			/**/, KmerDistanceCache & cachedCalculator
#endif
			) :
			KmerDistanceCalculator(
			/**/	dist,
			/**/	matrixId,
			/**/	kmerMode,
			/**/	fragMode,
			/**/	alphabet,
			/**/	fragLength,
			/**/	interval,
			/**/	kmerLength
			) {
			function<void(vector<Ranking>& rankings)> defaultQueryComplete = [&](vector<Ranking>& rankings) {};
			function<void(void)> defaultRunComplete = [&]() {};

			SetQueryComplete(defaultQueryComplete);
			SetRunComplete(defaultRunComplete);
		}

		virtual ~KmerSequenceRanker() {}

		double Threshold() {
			return threshold;
		}

		void SetThreshold(double value) {
			threshold = value;
		}

		function<void(vector<Ranking>& rankings)> & QueryComplete(void) {
			return queryComplete;
		}

		void SetQueryComplete(function<void(vector<Ranking>& rankings)> & value) {
			queryComplete = value;
		}

		function<void(void)> RunComplete(void) {
			return runComplete;
		}

		void SetRunComplete(function<void(void)> value) {
			runComplete = value;
		}

		void setup(void) {}

		int queryIdx = 0;

		void preProcess(void) {
			start_time =  time(NULL);
			Console::Error("Starting query " + Convert<size_t>::ToString(queryIdx + 1) + " at " + Ulong::ToString((ulong) start_time));
			rankings.clear();
			rankings.resize(subjectDataset->size());
		}

		void process(FastaSequence * querySeq, EncodedFastaSequence * subjectSeq, int subjectIdx, double distance) {
			Ranking ranking( querySeq->IdStr(), subjectSeq->IdStr(), distance, 0, 0);
			ranking.SetQueryId(querySeq->IdStr());
			ranking.SetQueryClass(querySeq->ClassLabel());
			rankings[subjectIdx] = ranking;
		}

		time_t start_time;

#if defined(USE_AMP_PROJECTOR)

		void postProcess( vector<int> & rowTotals, vector<int> & colTotals ) {
			{
				time_t end_time = time(NULL);

				Console::Error("Ending   query " 
					+ Convert<size_t>::ToString(queryIdx + 1) 
					+ " at " 
					+ Ulong::ToString((ulong) end_time)
					+ "("
					+ Double::ToString( (double) rankings.size() / ( end_time - start_time ) )
					+ " cmp/sec)"
					);
			}
			size_t queryKmerCount = (*queryDataset)[queryIdx]->Sequence().size() - kmerLength + 1;
			queryIdx++;

			for ( int i = 0; i < rankings.size(); i++ ) {
				size_t subjectKmerCount = (*subjectDataset)[i]->Sequence().size() - kmerLength + 1;
				auto distance = ((double) rowTotals[i] / queryKmerCount + (double) colTotals[i] / subjectKmerCount) / 2;
				rankings[i].SetDistance( distance );
			}

			sort(rankings.begin(), rankings.end(), Ranking::AscendingDistance);

			for ( size_t i = 0, len = rankings.size(); i < len; i++ ) {
				rankings[i].SetRank((int) i + 1);
			}

			queryComplete(rankings);
		}
#else
		void postProcess(void) {
			{
				time_t end_time = time(NULL);

				Console::Error("Ending   query "
					+ Convert<size_t>::ToString(queryIdx + 1)
					+ " at "
					+ Ulong::ToString((ulong) end_time)
					+ "("
					+ Double::ToString((double) rankings.size() / (end_time - start_time))
					+ " cmp/sec)"
					);
			}
			queryIdx++;

			sort(rankings.begin(), rankings.end(), Ranking::AscendingDistance);

			for ( size_t i = 0, len = rankings.size(); i < len; i++ ) {
				rankings[i].SetRank((int) i + 1);
			}

			queryComplete(rankings);
		}
#endif
		void cleanup(void) {
			runComplete();
		}


	};
}
