#pragma once

/**
**	This is the Dynamic Programming version of Similarity Projection.
**
**	The program exhaustively compares each query with the entire database
**	one sequence at a time, and uses exact kmer distances, subject to the 
**	choice of encoding and similarity matrix.
**
**	I expect it to be slow and accurate (competing with Smith-Waterman).
**	It is also able to do the bewildering combination for fragmentation
**	combination modes which are desirable for the first part of the thesis, 
**	where I explore the 
**	
**	Motivation for separating it is to allow me to cleanly implement the 
**	fission into (similarity lookup) vs (dual embedding) vs (dynamic programming) vs 
**	(VQ search) vs (Issl search).
*/

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <ctime>
#include <cfloat>
#include <omp.h>
#include <numeric>
#include <unordered_set>
#include <cstdint>

#include "Alphabet.hpp"
#include "Console.hpp"
#include "Delegates.hpp"
#include "FastaSequence.hpp"
#include "FragDistProcessor.hpp"
#include "FreeList.hpp"
#include "KmerCodebook.hpp"
#include "KmerSequenceRanker_Params.hpp"
#include "HausdorffCalculator.hpp"
#include "BestHitDetector.hpp"
#include "Projector.hpp"
#include "ProjectorBitEmbedding.hpp"
#include "ProjectorSlice.hpp"
#include "Util.hpp"
#include "Ranking.hpp"
#include "kNearestNeighbours.hpp"

#include "OmpTimer.h"

#if defined(USE_AMP)
#include "ProjectorAMP.hpp"
#endif

#if defined(USE_BIG_CACHE)
#include "KmerDistanceCache.hpp"
#endif


namespace QutBio {

	template <typename DistanceFunction, typename KmerType>
	class KmerDistanceCalculatorDP {
	protected:
		using D = DistanceFunction;
		using K = KmerType;
		using Cluster = KmerCluster<D, K>;

		/*
		**	Address of the current collection of pending rankings.
		**	(*rankings)[t] is the vector of rankings discovered by thread
		**	t.
		*/
		Rankings * rankings;

		/*
		**	Number of threads
		*/
		int numThreads;
		function<void(vector<Ranking> & rankings)> queryComplete;
		function<void(void)> runComplete;

		Alphabet * alphabet;
		SimilarityMatrix * matrix;
		FragmentAggregationMode * kmerMode;
		FragmentAggregationMode * fragMode;

		// Values needed by the Hausdorff family of distance calculators.
		uint fragmentLength;
		uint interval;
		uint kmerLength;
		const vector<EncodedFastaSequence *> * queryDataset;
		const vector<EncodedFastaSequence *> * subjectDataset;
		bool pushKmerDistances;
		vector<FreeList<SequenceDistanceFunction>> distanceFunctions;
		size_t maxQueryKmerCount;
		size_t maxQueryFragCount;

		DistanceFunction & distanceFunction;

		int skip;
		int maxRecords;

#if defined(WANT_DIAGNOSTIC_STREAM)
		ostream *diagnosticStream = 0;
#endif

	public:
		KmerDistanceCalculatorDP(
			SimilarityMatrix * matrix,
			FragmentAggregationMode * kmerMode,
			FragmentAggregationMode * fragMode,
			Alphabet * alphabet,
			uint fragLength,
			uint kmerLength,
			bool pushKmerDistances,
			DistanceFunction & distanceFunction,
			int skip,
			int maxRecords
		) :
			alphabet(alphabet),
			matrix(matrix),
			kmerMode(kmerMode),
			fragMode(fragMode),
			fragmentLength(fragLength),
			interval(1),
			kmerLength(kmerLength),
			pushKmerDistances(pushKmerDistances),
			distanceFunction(distanceFunction),
			skip(skip),
			maxRecords(maxRecords) {

#pragma omp parallel
			numThreads = omp_get_num_threads();

			distanceFunctions.resize(numThreads);
			pendingRankingsPerQuery.resize(numThreads);

			function<void(vector<Ranking>& rankings)> defaultQueryComplete = [&](vector<Ranking>& rankings) {};
			function<void(void)> defaultRunComplete = [&]() {};

			SetQueryComplete(defaultQueryComplete);
			SetRunComplete(defaultRunComplete);
		}

		virtual ~KmerDistanceCalculatorDP() {
			if (calculator) {
				delete calculator;
			}
		}

		function<void(vector<Ranking>& rankings)> & QueryComplete(void) {
			return queryComplete;
		}

		void SetQueryComplete(function<void(vector<Ranking> & rankings)> & value) {
			queryComplete = value;
		}

		function<void(void)> RunComplete(void) {
			return runComplete;
		}

		void SetRunComplete(function<void(void)> value) {
			runComplete = value;
		}

		/// <summary>
		/// Called once before the start of all processing.
		/// </summary>

		virtual void Setup(void) = 0;

		/// <summary>
		/// Called before each query is processed against the database.
		/// </summary>

		virtual void PreProcess(EncodedFastaSequence & querySeq, int queryIdx) = 0;

		/// <summary>
		/// Called to process each (query seq, db seq) pair.
		/// </summary>

		virtual void Process( EncodedFastaSequence & querySeq, size_t queryIdx, EncodedFastaSequence & subjectSeq, size_t subjectIdx, double distance) = 0;

		/// <summary>
		/// Called after all db sequences have been processed for given query sequence.
		/// </summary>

		virtual void PostProcess( EncodedFastaSequence & querySeq, int queryIdx) = 0;

		/// <summary>
		/// Called once after all queries have been processed.
		/// </summary>

		virtual void Cleanup(void) = 0;

		typedef KnnHeap<Ranking> KnnRankings;

		/// Run the job.
		void RunJob(
			vector<EncodedFastaSequence *> & query,
			vector<EncodedFastaSequence *> & db
		) {
			queryDataset = &query;
			subjectDataset = &db;

			int maxQueryLength = numeric_limits<int>::min();

			for (auto seq : query) {
				int len = (int)seq->Sequence().size();

				if (len > maxQueryLength) maxQueryLength = len;
			}

			int maxSubjectLength = numeric_limits<int>::min();

			for (auto seq : db) {
				int len = (int)seq->Sequence().size();

				if (len > maxSubjectLength) maxSubjectLength = len;
			}

			maxQueryKmerCount = maxQueryLength - kmerLength + 1;
			maxQueryFragCount = Fragment::GetCount(maxQueryKmerCount, fragmentLength);

			Setup();

			// Later we intend to sort the rankings, but sort is single-threaded.
			// So... accumulate rankings until we have enough to give a copy to each 
			// thread for sorting.
			vector<Rankings> pendingRankings(numThreads);

			for (auto & rankings : pendingRankings) {
				rankings.resize(numThreads);
			}

			int numPending = 0;

			KnnRankings knnAll(maxRecords);
			vector<KnnRankings> knnPerThread(numThreads);

			for (int i = 0; i < numThreads; i++) {
				knnPerThread[i].setCapacity(maxRecords);
			}


			for (size_t queryIdx = 0; queryIdx < query.size(); queryIdx++) {
				rankings = &pendingRankings[numPending++];

				auto & querySeq = *query[queryIdx];
				PreProcess(querySeq, (int)queryIdx);
				DoPairwise(queryIdx, querySeq, db, maxQueryLength, maxSubjectLength, knnPerThread, knnAll);
				PostProcess(querySeq, (int)queryIdx);

				if (numPending == numThreads) {
					BatchProcessRankings(pendingRankings, numPending);
					numPending = 0;
				}
			}

			if (numPending > 0) {
				BatchProcessRankings(pendingRankings, numPending);
			}

			Cleanup();
		}

		vector<vector<Ranking>> pendingRankingsPerQuery;

		/**
		 **	Processes pending rankings collections a set of numThreads queries.
		 **	The sort is done in parallel, followed by a serial write operation.
		 **/

		void BatchProcessRankings(
			vector<Rankings> &pendingRankings,
			int numberRemaining
		) {
#pragma omp parallel for ordered 
			for (int i = 0; i < numberRemaining; i++) {
				Rankings &rankingsPerThread = pendingRankings[i];
				vector<Ranking> &rankings = pendingRankingsPerQuery[i];
				rankings.clear();

				for (int t = 0; t < numThreads; t++) {
					rankings.insert(rankings.end(), rankingsPerThread[t].begin(), rankingsPerThread[t].end());
				}

				sort(rankings.begin(), rankings.end(), Ranking::AscendingDistance);

#pragma omp ordered
				queryComplete(pendingRankingsPerQuery[i]);
			}
		}

		const SimilarityMatrix & Matrix() {
			return *matrix;
		}

		vector<vector<uint16_t>> colMinima;
		vector<vector<uint16_t>> rowMinima;
		vector<FragDistProcessor_Frag<SubVector<uint16_t>, uint16_t>> fragDistProcessor_frag;
		vector<FragDistProcessor_NoFrag<SubVector<uint16_t>, uint16_t>> fragDistProcessor_noFrag;
		vector<bool> subjectSeen;
		vector<vector<size_t>> subjectOffset;
		vector<size_t> queryFragMapping;
		vector<size_t> rankingOffsetPerThread;

		typedef pair<size_t, Cluster *> OffsetClusterMapping;

		vector<vector<OffsetClusterMapping>> relevantClustersPerThread;

		template<typename vectorType>
		double FragHausdorffAverageAverage_FragDistCache(
			vectorType & rowMinima,
			vectorType & colMinima,
			size_t queryFragCount,
			size_t queryKmerCount,
			size_t subjectFragCount,
			size_t subjectKmerCount,
			int fragmentLength,
			int shift,
			int querySkip
		) {
			int totXY = 0;
			int totYX = 0;

			int queryWeight = 0;

			// Compute average one-way distance from query to subject, weighting 
			//	each fragment by the number of kmers contained therein.
			if (querySkip > 1 && fragmentLength == 1) {
				for (size_t i = 0; i < queryKmerCount; i += querySkip) {
					totXY += rowMinima[i] << shift;
					queryWeight++;
				}
			}
			else {
				int queryQuotient = (int)queryKmerCount / fragmentLength;
				int queryRemainder = (int)queryKmerCount % fragmentLength;

				queryWeight = queryQuotient * fragmentLength + queryRemainder;

				// Compute average one-way distance from query to subject, weighting 
				//	each fragment by the number of kmers contained therein.
				for (int j = 0; j < queryQuotient; j++) {
					totXY += rowMinima[j] << shift;
				}

				totXY *= fragmentLength;

				if (queryRemainder) totXY += rowMinima[queryQuotient] * queryRemainder;
			}

			int subjectQuotient = (int)subjectKmerCount / fragmentLength;
			int subjectRemainder = (int)subjectKmerCount % fragmentLength;
			int subjectWeight = subjectQuotient * fragmentLength + subjectRemainder;

			// Compute average one-way distance from subject to query, weighting 
			//	each fragment by the number of kmers contained therein.
			for (int j = 0; j < subjectQuotient; j++) {
				totYX += colMinima[j] << shift;
			}

			totYX *= fragmentLength;

			if (subjectRemainder) totYX += colMinima[subjectQuotient] * subjectRemainder;

			double avgXY = (double)totXY / queryWeight;
			double avgYX = (double)totYX / subjectWeight;

			double result = (avgXY + avgYX) / 2;
			return result;
		}

		void DoPairwise(
			size_t queryIdx,
			EncodedFastaSequence &querySeq,
			vector<EncodedFastaSequence *> & db,
			size_t maxQueryLength,
			size_t maxSubjectLength,
			vector<KnnRankings> &knnPerThread,
			KnnRankings &knnAll
		) {
			for (KnnRankings & knn : knnPerThread) {
				knn.clear();
			}

#pragma omp parallel for
			for (int64_t subjectIdx = 0; subjectIdx < (int64_t)db.size(); subjectIdx++) {
				int thread = omp_get_thread_num();
				KnnRankings & knn = knnPerThread[thread];
				FreeList<SequenceDistanceFunction> & sequenceDistances = distanceFunctions[thread];
				EncodedFastaSequence & subjectSeq = *db[subjectIdx];
				DoPairwise(queryIdx, querySeq, maxQueryLength, subjectIdx, subjectSeq, maxSubjectLength, knn, sequenceDistances);
			}

			// Now get the maxRecords nearest overall.
			{
				knnAll.clear();

				for (KnnRankings & knn : knnPerThread) {
					for (Ranking & r : knn) {
						knnAll.push(r);
					}
				}
			}

			// Copy maxRecords nearest to the pending rankings vector.
			// This stuff really needs to move to another program!
			{
				// This function uses the Knn objects to capture a small subset.
				// so we only end up with a single thread to copy rankings data. 
				for (size_t i = 0; i < rankings->size(); i++) {
					(*rankings)[i].clear();
				}

				(*rankings)[0].insert((*rankings)[0].end(), knnAll.begin(), knnAll.end());
			}
		}

		HausdorffCalculator *calculator = 0;

		void DoPairwise(
			size_t queryIdx,
			EncodedFastaSequence & querySeq,
			size_t maxQueryLength,
			size_t subjectIdx,
			EncodedFastaSequence & subjectSeq,
			size_t maxSubjectLength,
			KnnHeap<Ranking> & knn,
			FreeList<SequenceDistanceFunction> & distanceFunctions
		) {
			if (!calculator) {
				calculator = new HausdorffCalculator(
					matrix,
					kmerLength,
					kmerMode,
					fragMode,
					alphabet,
					fragmentLength,
					maxQueryLength,
					maxSubjectLength
				);
			}

			double distance = calculator->ComputeDistance(*querySeq.base, *subjectSeq.base);

			Ranking r(querySeq.IdStr(), subjectSeq.IdStr(), distance, 0, 0);
			knn.push(r);

#if defined(WANT_DIAGNOSTIC_STREAM) && WANT_DIAGNOSTIC_STREAM
			if (diagnosticStream) {
#pragma omp critical
				{
					ostream & out(*diagnosticStream);

					out << querySeq.IdStr() << "," << subjectSeq.IdStr() << "," << distance << "," << 0 << "\n";

					out << querySeq.IdStr();

					for (size_t i = 0; i < calculator->queryFragCount; i++) {
						int distance = calculator->rowMinima[i];
						out << "," << distance;
					}

					out << "\n" << subjectSeq.IdStr();

					for (size_t i = 0; i < calculator->subjectFragCount; i++) {
						int distance = calculator->colMinima[i];
						out << "," << distance;
					}

					out << "\n" << querySeq.IdStr() << "," << subjectSeq.IdStr();

					for (size_t i = 0; i < calculator->queryFragCount; i++) {
						out << ";";

						for (size_t j = 0; j < calculator->subjectFragCount; j++) {
							if (j > 0) {
								out << ",";
							}
							out << calculator->kmerDistCache(i,j);
						}
					}

					out << "\n";
				}
			}
#endif
		}
	};
}
