#pragma once

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
#include "Projector.hpp"
#include "ProjectorBitEmbedding.hpp"
#include "ProjectorSlice.hpp"
#include "Util.hpp"
#include "Ranking.hpp"
#include "kNearestNeighbours.hpp"
#include "Exception.hpp"

#include "OmpTimer.h"

#if defined(USE_AMP)
#include "ProjectorAMP.hpp"
#endif

#if defined(USE_BIG_CACHE)
#include "KmerDistanceCache.hpp"
#endif


namespace QutBio {

	template <typename DistanceFunction, typename KmerType>
	class KmerDistanceCalculator {
	protected:
		using D = DistanceFunction;
		using K = KmerType;
		using Cluster = KmerCluster<D, K>;
		using Codebook = KmerCodebook<D, K>;

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
		function<void( vector<Ranking> & rankings )> queryComplete;
		function<void( void )> runComplete;

		Alphabet * alphabet;
		SimilarityMatrix * matrix;
		FragmentAggregationMode * kmerMode;
		FragmentAggregationMode * fragMode;

		// Values needed by the Hausdorff family of distance calculators.
		uint fragmentLength;
		uint interval;
		uint kmerLength;
		vector<EncodedFastaSequence *> * queryDataset;
		vector<EncodedFastaSequence *> * subjectDataset;
		Distance thresholdDistance = BAD_DIST;
		Distance defaultDistance = BAD_DIST;
		Codebook * codebook = 0;
		bool pushKmerDistances;
		vector<FreeList<SequenceDistanceFunction>> distanceFunctions;
		size_t maxQueryKmerCount;
		size_t maxQueryFragCount;

		DistanceFunction & distanceFunction;

		int skip;
		int maxRecords;

#if WANT_DIAGNOSTIC_STREAM
		ostream *diagnosticStream = 0;
#endif

	public:
		KmerDistanceCalculator(
			SimilarityMatrix * matrix,
			FragmentAggregationMode * kmerMode,
			FragmentAggregationMode * fragMode,
			Alphabet * alphabet,
			uint fragLength,
			uint kmerLength,
			Distance thresholdDistance,
			Distance defaultDistance,
			bool pushKmerDistances,
			DistanceFunction & distanceFunction,
			int skip,
			int maxRecords
		) :
			alphabet( alphabet ),
			matrix( matrix ),
			kmerMode( kmerMode ),
			fragMode( fragMode ),
			fragmentLength( fragLength ),
			interval( 1 ),
			kmerLength( kmerLength ),
			thresholdDistance( thresholdDistance ),
			defaultDistance( defaultDistance ),
			pushKmerDistances( pushKmerDistances ),
			distanceFunction( distanceFunction ),
			skip( skip ),
			maxRecords( maxRecords ) {
#if USE_OMP
#pragma omp parallel
#endif
			numThreads = omp_get_num_threads();

			distanceFunctions.resize( numThreads );
			pendingRankingsPerQuery.resize( numThreads );

			function<void( vector<Ranking>& rankings )> defaultQueryComplete = [&]( vector<Ranking>& rankings ) {};
			function<void( void )> defaultRunComplete = [&]() {};

			SetQueryComplete( defaultQueryComplete );
			SetRunComplete( defaultRunComplete );
		}

		virtual ~KmerDistanceCalculator() {}

		void SetThreshold( Distance thresholdDistance, Distance defaultDistance ) {
			assert_false( IS_BAD_DIST( thresholdDistance ) );
			assert_false( IS_BAD_DIST( defaultDistance ) );
			assert_true( thresholdDistance <= defaultDistance );

			thresholdDistance = thresholdDistance;
			defaultDistance = defaultDistance;
		}

		double ThresholdDistance() { return KmerDistanceCalculator<D, K>::thresholdDistance; }

		double DefaultDistance() { return KmerDistanceCalculator<D, K>::defaultDistance; }

		function<void( vector<Ranking>& rankings )> & QueryComplete( void ) {
			return queryComplete;
		}

		void SetQueryComplete( function<void( vector<Ranking> & rankings )> & value ) {
			queryComplete = value;
		}

		function<void( void )> RunComplete( void ) {
			return runComplete;
		}

		void SetRunComplete( function<void( void )> value ) {
			runComplete = value;
		}

		/// <summary>
		/// Called once before the start of all processing.
		/// </summary>

		virtual void setup( void ) = 0;

		/// <summary>
		/// Called before each query is processed against the database.
		/// </summary>

		virtual void preProcess( EncodedFastaSequence & querySeq, int queryIdx ) = 0;

		/// <summary>
		/// Called to process each (query seq, db seq) pair.
		/// </summary>

		virtual void process( EncodedFastaSequence & querySeq, size_t queryIdx, EncodedFastaSequence & subjectSeq, size_t subjectIdx, double distance ) = 0;

		/// <summary>
		/// Called after all db sequences have been processed for given query sequence.
		/// </summary>

		virtual void postProcess( EncodedFastaSequence & querySeq, int queryIdx ) = 0;

		/// <summary>
		/// Called once after all queries have been processed.
		/// </summary>

		virtual void cleanup( void ) = 0;

		typedef KnnHeap<Ranking> KnnRankings;

		/// Run the job.
		void RunJob(
			vector<EncodedFastaSequence *> & query,
			vector<EncodedFastaSequence *> & db
		) {
			queryDataset = &query;
			subjectDataset = &db;

			int maxQueryLength = numeric_limits<int>::min();

			for ( auto seq_ : query ) {
				auto & seq = *seq_;
				int len = (int) seq.Sequence().size();

				if ( len > maxQueryLength ) maxQueryLength = len;
			}

			int maxSubjectLength = numeric_limits<int>::min();

			for ( auto seq_ : db ) {
				auto & seq = *seq_;
				int len = (int) seq.Sequence().size();

				if ( len > maxSubjectLength ) maxSubjectLength = len;
			}

			maxQueryKmerCount = maxQueryLength - kmerLength + 1;
			maxQueryFragCount = Fragment::GetCount( maxQueryKmerCount, fragmentLength );

			setup();

			// Later we intend to sort the rankings, but sort is single-threaded.
			// So... accumulate rankings until we have enough to give a copy to each 
			// thread for sorting.
			vector<Rankings> pendingRankings( numThreads );

			for ( auto & rankings : pendingRankings ) {
				rankings.resize( numThreads );
			}

			int numPending = 0;

			KnnRankings knnAll( maxRecords );

			vector<KnnRankings> knnPerThread( numThreads );

			for ( int i = 0; i < numThreads; i++ ) {
				knnPerThread[i].setCapacity( maxRecords );
			}

			vector<HausdorffCalculator *> calcPerThread( numThreads );

			for ( int i = 0; i < numThreads; i++ ) {
				calcPerThread[i] = new HausdorffCalculator(
					matrix,
					kmerLength,
					kmerMode,
					fragMode,
					alphabet,
					fragmentLength,
					maxQueryLength,
					maxSubjectLength
				);

				if ( !IS_BAD_DIST( thresholdDistance ) ) {
					calcPerThread[i]->SetThreshold( thresholdDistance, defaultDistance );
				}
			}

			for ( size_t queryIdx = 0; queryIdx < query.size(); queryIdx++ ) {
				rankings = &pendingRankings[numPending++];

				EncodedFastaSequence & querySeq = *query[queryIdx];
				preProcess( querySeq, (int) queryIdx );

				if ( codebook && !IS_BAD_DIST( thresholdDistance ) ) {
					PopulateFragDistCache( querySeq, queryIdx, db );
				}
				else {
					DoPairwise( queryIdx, querySeq, db, maxQueryLength, maxSubjectLength, knnPerThread, calcPerThread, knnAll );
				}

				postProcess( querySeq, (int) queryIdx );

				if ( numPending == numThreads ) {
					batchProcessRankings( pendingRankings, numPending );
					numPending = 0;
				}
			}

			for ( int i = 0; i < numThreads; i++ ) {
				delete calcPerThread[i];
			}

			if ( numPending > 0 ) {
				batchProcessRankings( pendingRankings, numPending );
			}

			cleanup();
		}

		vector<vector<Ranking>> pendingRankingsPerQuery;

		/**
		**	Processes pending rankings collections a set of numThreads queries.
		**	The sort is done in parallel, followed by a serial write operation.
		**/

		void batchProcessRankings(
			vector<Rankings> &pendingRankings,
			int numberRemaining
		) {
#if USE_OMP
#pragma omp parallel for ordered 
#endif
			for ( int i = 0; i < numberRemaining; i++ ) {

				Rankings &rankingsPerThread = pendingRankings[i];
				vector<Ranking> &rankings = pendingRankingsPerQuery[i];
				rankings.clear();

				for ( int t = 0; t < numThreads; t++ ) {
					rankings.insert( rankings.end(), rankingsPerThread[t].begin(), rankingsPerThread[t].end() );
				}

				sort( rankings.begin(), rankings.end(), Ranking::AscendingDistance );

#if USE_OMP
#pragma omp ordered
#endif
				queryComplete( pendingRankingsPerQuery[i] );
			}
		}

		const SimilarityMatrix & Matrix() {
			return *matrix;
		}

		using VectorType = vector<Distance>;

		vector<VectorType> colMinima;
		vector<VectorType> rowMinima;
		vector<FragDistProcessor_Frag  <MutableSubVector<Distance>, Distance>> fragDistProcessor_frag;
		vector<FragDistProcessor_NoFrag<MutableSubVector<Distance>, Distance>> fragDistProcessor_noFrag;
		vector<bool> subjectSeen;
		vector<vector<size_t>> subjectOffset;
		vector<size_t> queryFragMapping;
		vector<size_t> rankingOffsetPerThread;

		typedef pair<size_t, Cluster *> OffsetClusterMapping;

		// TODO: Replace this with a collection of priority queues
		vector<vector<OffsetClusterMapping>> relevantClustersPerThread;

		void PopulateFragDistCache(
			EncodedFastaSequence & querySeq,
			size_t queryIdx,
			vector<EncodedFastaSequence *> & db
		) {
			// TODO: move the isFirst stuff into the constructor.
			// TODO: split this into two classes because the pairwise version needs different parallel infrastructure. 
			static bool auxDataDone = false;
			bool isFirst = false;

			// ( cerr << "A" << endl ).flush();

			if ( !auxDataDone ) {
				auxDataDone = true;
				isFirst = true;
			}

			const size_t dbSize = db.size();
			const int T = numThreads;
			const int N = (int) codebook->Size();

			// TODO: Convert Shift distance into a parameter or get rid of it.
			const int shift_distance = 0;

			const size_t maxPartSize = (dbSize + numThreads - 1) / numThreads;

			if ( isFirst ) {
				// Move this into an initialisation function elsewhere
				fragDistProcessor_frag.resize( dbSize );
				colMinima.resize( numThreads );
				rowMinima.resize( numThreads );
				subjectOffset.resize( numThreads );
				vector<size_t> totLength( numThreads, 0 );
				queryFragMapping.resize( maxQueryKmerCount );
				rankingOffsetPerThread.resize( numThreads );
#if USE_OMP
#pragma omp parallel
#endif
				{
					int t = omp_get_thread_num();
					subjectOffset[t].resize( maxPartSize );
					rankings[t].reserve( maxRecords );

					for ( size_t n = 0; n < maxPartSize; n++ ) {
						size_t i = t * maxPartSize + n;

						if ( i >= dbSize ) break;

						auto & seq = *db[i];
						auto & subjectStr = seq.Sequence();
						size_t subjectKmerCount = subjectStr.size() - kmerLength + 1;
						size_t subjectFragCount = Fragment::GetCount( subjectKmerCount, fragmentLength );

						subjectOffset[t][n] = totLength[t];
						totLength[t] += subjectFragCount;
					}

					colMinima[t].resize( totLength[t] );
					rowMinima[t].resize( maxPartSize * maxQueryFragCount );

					for ( size_t n = 0; n < maxPartSize; n++ ) {
						size_t i = t * maxPartSize + n;

						if ( i >= dbSize ) break;

						auto & seq = *db[i];
						auto & subjectStr = seq.Sequence();
						auto subjectKmerCount = subjectStr.size() - kmerLength + 1;
						auto subjectFragCount = Fragment::GetCount( subjectKmerCount, fragmentLength );

						MutableSubVector<Distance> thisColMinima( &colMinima[t], subjectOffset[t][n], subjectFragCount );
						MutableSubVector<Distance> thisRowMinima( &rowMinima[t], n * maxQueryFragCount, maxQueryFragCount );

						auto & f( fragDistProcessor_frag[i] );
						f.rowMinima = thisRowMinima;
						f.colMinima = thisColMinima;
						f.queryKmerCount = maxQueryKmerCount;
						f.subjectKmerCount = subjectKmerCount;
						f.queryFragCount = maxQueryFragCount;
						f.subjectFragCount = subjectFragCount;
						f.fragmentLength = fragmentLength;
						f.fragsPerTile = 1; // fragsPerTile;
					}
				}

				subjectSeen.resize( db.size() );

				relevantClustersPerThread.resize( T );

				for ( int t = 0; t < T; t++ ) {
					relevantClustersPerThread[t].reserve( maxPartSize );
				}
			}

			size_t queryKmerCount = querySeq.Sequence().size() - kmerLength + 1;
			size_t queryFragCount = Fragment::GetCount( queryKmerCount, fragmentLength );

			for ( int i = 0; i < numThreads; i++ ) {
				relevantClustersPerThread[i].clear();
				(*rankings)[i].clear();
				rankingOffsetPerThread[i] = 0;
			}

#if USE_OMP
#pragma omp parallel for
#endif
			for ( size_t i = 0; i < db.size(); i++ ) {
				fragDistProcessor_frag[i].queryFragCount = queryFragCount;
				fragDistProcessor_frag[i].queryKmerCount = queryKmerCount;
				subjectSeen[i] = false;
			}

#if USE_OMP
#pragma omp parallel
#endif
			{
				int t = omp_get_thread_num();
				int start = (N * t + T / 2) / T;
				int end = (N * (t + 1) + T / 2) / T;

				for ( int queryPos = 0; queryPos < (int) queryKmerCount; queryPos += skip ) {
					auto queryWord = querySeq.GetEncodedKmer( queryPos );

					for ( int codebookPos = start; codebookPos < end; codebookPos++ ) {
						EncodedKmer codeword = (EncodedKmer) codebook->kmerData.row( codebookPos );

						if ( distanceFunction( queryWord, codeword, kmerLength ) <= thresholdDistance ) {
							relevantClustersPerThread[t].emplace_back( queryPos, codebook->codebook[codebookPos] );
						}
					}
				}
			}

#if WANT_DIAGNOSTIC_STREAM
			if ( diagnosticStream ) {
#if USE_OMP
#pragma omp critical
#endif
				{
					ostream & str = (*diagnosticStream);
					str << querySeq.IdStr() << "," << "Selected Clusters\n";

					for ( int t = 0; t < numThreads; t++ ) {
						for ( OffsetClusterMapping & tuple : relevantClustersPerThread[t] ) {
							size_t queryPos = tuple.first;
							Cluster * cluster = tuple.second;
							auto & protoInstance = cluster->prototype->Instances()[0];
							str << querySeq.IdStr() << "," << queryPos << "," << protoInstance.sequence.IdStr() << "," << protoInstance.kmerPosition << "\n";
						}
					}

					str << "\n\n";
				}
			}
#endif

#if USE_OMP
#pragma omp parallel
#endif
			{
				int threadId = omp_get_thread_num();

				for ( int t = 0; t < numThreads; t++ ) {
					for ( auto & hit : relevantClustersPerThread[t] ) {
						Cluster & cluster = *hit.second;
						auto & kmers = cluster.kmersPerThread[threadId];

						size_t queryPos = hit.first;
						auto queryWord = querySeq.GetEncodedKmer( queryPos );
						Distance distance = numeric_limits<Distance>::max();

						for ( auto & kmer : kmers ) {
							KmerWord * subjectWord = kmer.first->PackedEncoding();
							distance = distanceFunction( queryWord, subjectWord, kmerLength );

							if ( distance < defaultDistance ) {
								for ( auto & instance : kmer.first->Instances() ) {
									auto & subjectSeq = instance.sequence;
									auto subjectIdx = subjectSeq.position;
									auto & fragDistProcessor = fragDistProcessor_frag[subjectIdx];

									if ( !subjectSeen[subjectIdx] ) {
										fragDistProcessor.Reset( queryFragCount, defaultDistance );
										subjectSeen[subjectIdx] = true;
									}

									fragDistProcessor.ProcessDistance( queryPos, instance.kmerPosition, distance );
								}
							}
						}
					}
				}

			}

			if ( fragMode == FragmentAggregationMode::HausdorffAverageAverage() ) {
				HausdorffProcess( db, queryFragCount, queryKmerCount, shift_distance, querySeq, FragHausdorffAverageAverage_FragDistCache );
			}
			else if ( fragMode == FragmentAggregationMode::BestOfBest() ) {
				HausdorffProcess( db, queryFragCount, queryKmerCount, shift_distance, querySeq, FragBestOfBest_FragDistCache );
			}
			else if ( fragMode == FragmentAggregationMode::HausdorffAverage() ) {
				HausdorffProcess( db, queryFragCount, queryKmerCount, shift_distance, querySeq, FragHausdorffAverage_FragDistCache );
			}
			else if ( fragMode == FragmentAggregationMode::Hausdorff() ) {
				HausdorffProcess( db, queryFragCount, queryKmerCount, shift_distance, querySeq, FragHausdorff_FragDistCache );
			}
			else {
				HausdorffProcess( db, queryFragCount, queryKmerCount, shift_distance, querySeq, FragHausdorffAverageAverage_FragDistCache );
			}
		}

		using Processor = function<double (
			MutableSubVector<Distance> & rowMinima,
			MutableSubVector<Distance> & colMinima,
			size_t queryFragCount,
			size_t queryKmerCount,
			size_t subjectFragCount,
			size_t subjectKmerCount,
			int fragmentLength,
			int shift,
			int querySkip
		)>;

		void HausdorffProcess(
			std::vector<QutBio::EncodedFastaSequence *> & db,
			size_t &queryFragCount,
			size_t &queryKmerCount,
			const int &shift_distance,
			QutBio::EncodedFastaSequence & querySeq,
			Processor process
		) {
#if USE_OMP
#pragma omp parallel for
#endif
			for ( int64_t subjectIdx = 0; subjectIdx < (int64_t) db.size(); subjectIdx++ ) {
				if ( !subjectSeen[subjectIdx] ) continue;

				auto & subjectSeq = *db[subjectIdx];
				int t = omp_get_thread_num();

				auto & fragDistProcessor = fragDistProcessor_frag[subjectIdx];

				// TODO: Move this into the FragDistProcessor, where I'm sure there is enough information to handle it.
				double distance = process(
					fragDistProcessor.rowMinima,
					fragDistProcessor.colMinima,
					queryFragCount,
					queryKmerCount,
					fragDistProcessor.subjectFragCount,
					fragDistProcessor.subjectKmerCount,
					fragmentLength,
					shift_distance,
					skip
				);

				(*rankings)[t].emplace_back( querySeq.IdStr(), subjectSeq.IdStr(), distance, 0, fragDistProcessor.hits );

#if WANT_DIAGNOSTIC_STREAM 
				if ( diagnosticStream ) {
#if USE_OMP
#pragma omp critical
#endif
					{
						ostream & str = (*diagnosticStream);
						str << querySeq.IdStr() << "," << subjectSeq.IdStr() << "," << distance << "," << fragDistProcessor.hits << "\n\n";

						str << querySeq.IdStr();

						for ( size_t i = 0; i < queryFragCount; i++ ) {
							str << "," << fragDistProcessor.rowMinima[i];
						}

						str << "\n\n" << subjectSeq.IdStr();

						for ( size_t i = 0; i < fragDistProcessor.subjectFragCount; i++ ) {
							str << "," << fragDistProcessor.colMinima[i];
						}

						str << "\n\n";
					}
				}
#endif
			}
		}

		static double FragHausdorffAverageAverage_FragDistCache(
			MutableSubVector<Distance> & rowMinima,
			MutableSubVector<Distance> & colMinima,
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
			if ( querySkip > 1 && fragmentLength == 1 ) {
				for ( size_t i = 0; i < queryKmerCount; i += querySkip ) {
					totXY += rowMinima[i] << shift;
					queryWeight++;
				}
			}
			else {
				int queryQuotient = (int) queryKmerCount / fragmentLength;
				int queryRemainder = (int) queryKmerCount % fragmentLength;

				queryWeight = queryQuotient * fragmentLength + queryRemainder;

				// Compute average one-way distance from query to subject, weighting 
				//	each fragment by the number of kmers contained therein.
				for ( int j = 0; j < queryQuotient; j++ ) {
					totXY += rowMinima[j] << shift;
				}

				totXY *= fragmentLength;

				if ( queryRemainder ) totXY += (rowMinima[queryQuotient] << shift) * queryRemainder;
			}

			int subjectQuotient = (int) subjectKmerCount / fragmentLength;
			int subjectRemainder = (int) subjectKmerCount % fragmentLength;
			int subjectWeight = subjectQuotient * fragmentLength + subjectRemainder;

			// Compute average one-way distance from subject to query, weighting 
			//	each fragment by the number of kmers contained therein.
			for ( int j = 0; j < subjectQuotient; j++ ) {
				totYX += colMinima[j] << shift;
			}

			totYX *= fragmentLength;

			if ( subjectRemainder ) totYX += (colMinima[subjectQuotient] << shift) * subjectRemainder;

			double avgXY = (double) totXY / queryWeight;
			double avgYX = (double) totYX / subjectWeight;

			double result = (avgXY + avgYX) / 2;
			return result;
		}

		static double FragHausdorffAverage_FragDistCache(
			MutableSubVector<Distance> & rowMinima,
			MutableSubVector<Distance> & colMinima,
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
			if ( querySkip > 1 && fragmentLength == 1 ) {
				for ( size_t i = 0; i < queryKmerCount; i += querySkip ) {
					totXY += rowMinima[i] << shift;
					queryWeight++;
				}
			}
			else {
				int queryQuotient = (int) queryKmerCount / fragmentLength;
				int queryRemainder = (int) queryKmerCount % fragmentLength;

				queryWeight = queryQuotient * fragmentLength + queryRemainder;

				// Compute average one-way distance from query to subject, weighting 
				//	each fragment by the number of kmers contained therein.
				for ( int j = 0; j < queryQuotient; j++ ) {
					totXY += rowMinima[j] << shift;
				}

				totXY *= fragmentLength;

				if ( queryRemainder ) totXY += (rowMinima[queryQuotient] << shift) * queryRemainder;
			}

			int subjectQuotient = (int) subjectKmerCount / fragmentLength;
			int subjectRemainder = (int) subjectKmerCount % fragmentLength;
			int subjectWeight = subjectQuotient * fragmentLength + subjectRemainder;

			// Compute average one-way distance from subject to query, weighting 
			//	each fragment by the number of kmers contained therein.
			for ( int j = 0; j < subjectQuotient; j++ ) {
				totYX += colMinima[j] << shift;
			}

			totYX *= fragmentLength;

			if ( subjectRemainder ) totYX += (colMinima[subjectQuotient] << shift) * subjectRemainder;

			double avgXY = (double) totXY / queryWeight;
			double avgYX = (double) totYX / subjectWeight;

			double result = std::max( avgXY, avgYX );
			return result;
		}

		static double FragHausdorff_FragDistCache(
			MutableSubVector<Distance> & rowMinima,
			MutableSubVector<Distance> & colMinima,
			size_t queryFragCount,
			size_t queryKmerCount,
			size_t subjectFragCount,
			size_t subjectKmerCount,
			int fragmentLength,
			int shift,
			int querySkip
		) {
			int maxXY = numeric_limits<int>::min();
			int maxYX = numeric_limits<int>::min();

			// Compute average one-way distance from query to subject, weighting 
			//	each fragment by the number of kmers contained therein.
			if ( querySkip > 1 && fragmentLength == 1 ) {
				for ( size_t i = 0; i < queryKmerCount; i += querySkip ) {
					int d = rowMinima[i] << shift;
					maxXY = std::max( maxXY, d );
				}
			}
			else {
				int queryQuotient = (int) queryKmerCount / fragmentLength;
				int queryRemainder = (int) queryKmerCount % fragmentLength;

				// Compute average one-way distance from query to subject.
				for ( int j = 0; j < queryQuotient; j++ ) {
					int d = rowMinima[j] << shift;
					maxXY = std::max( maxXY, d );
				}

				if ( queryRemainder ) {
					int d = rowMinima[queryQuotient] << shift;
					maxXY = std::max( maxXY, d );
				}
			}

			int subjectQuotient = (int) subjectKmerCount / fragmentLength;
			int subjectRemainder = (int) subjectKmerCount % fragmentLength;

			// Compute maximum one-way distance from subject to query.
			for ( int j = 0; j < subjectQuotient; j++ ) {
				int d = colMinima[j] << shift;
				maxYX = std::max( maxYX, d );
			}

			if ( subjectRemainder ) {
				int d = colMinima[subjectQuotient] << shift;
				maxYX = std::max( maxYX, d );
			}

			double result = std::max( maxXY, maxYX );
			return result;
		}

		static double FragBestOfBest_FragDistCache(
			MutableSubVector<Distance> & rowMinima,
			MutableSubVector<Distance> & colMinima,
			size_t queryFragCount,
			size_t queryKmerCount,
			size_t subjectFragCount,
			size_t subjectKmerCount,
			int fragmentLength,
			int shift,
			int querySkip
		) {
			auto result = rowMinima[0];

			if ( querySkip > 1 && fragmentLength == 1 ) {
				for ( size_t i = 0; i < queryKmerCount; i += querySkip ) {
					double d = rowMinima[i] << shift;
					if ( d < result ) result = d;
				}
			}
			else {
				int queryQuotient = (int) queryKmerCount / fragmentLength;
				int queryRemainder = (int) queryKmerCount % fragmentLength;

				for ( int j = 0; j < queryQuotient; j++ ) {
					double d = rowMinima[j] << shift;
					if ( d < result ) result = d;
				}

				if ( queryRemainder ) {
					double d = rowMinima[queryQuotient] << shift;
					if ( d < result ) result = d;
				}
			}

			int subjectQuotient = (int) subjectKmerCount / fragmentLength;
			int subjectRemainder = (int) subjectKmerCount % fragmentLength;

			// Compute average one-way distance from subject to query, weighting 
			//	each fragment by the number of kmers contained therein.
			for ( int j = 0; j < subjectQuotient; j++ ) {
				double d = colMinima[j] << shift;
				if ( d < result ) result = d;
			}

			if ( subjectRemainder ) {
				double d = colMinima[subjectQuotient] << shift;
				if ( d < result ) result = d;
			}

			return result;
		}

		void DoPairwise(
			size_t queryIdx,
			EncodedFastaSequence &querySeq,
			vector<EncodedFastaSequence *> & db,
			size_t maxQueryLength,
			size_t maxSubjectLength,
			vector<KnnRankings> &knnPerThread,
			vector<HausdorffCalculator *> calcPerThread,
			KnnRankings &knnAll
		) {
			for ( KnnRankings & knn : knnPerThread ) {
				knn.clear();
			}

#if USE_OMP
#pragma omp parallel for
#endif
			for ( int64_t subjectIdx = 0; subjectIdx < (int64_t) db.size(); subjectIdx++ ) {
				int thread = omp_get_thread_num();
				KnnRankings & knn = knnPerThread[thread];
				auto & subjectSeq = *db[subjectIdx];
				DoPairwise( queryIdx, querySeq, maxQueryLength, subjectIdx, subjectSeq, maxSubjectLength, knn, calcPerThread[thread] );
			}

			// Now get the maxRecords nearest overall.
			{
				knnAll.clear();

				for ( KnnRankings & knn : knnPerThread ) {
					for ( Ranking & r : knn ) {
						knnAll.push( r );
					}
				}
			}

			// Copy maxRecords nearest to the pending rankings vector.
			// This stuff really needs to move to another program!
			{
				// This function uses the Knn objects to capture a small subset.
				// so we only end up with a single thread to copy rankings data. 
				for ( size_t i = 0; i < rankings->size(); i++ ) {
					(*rankings)[i].clear();
				}

				(*rankings)[0].insert( (*rankings)[0].end(), knnAll.begin(), knnAll.end() );
			}
		}

		void DoPairwise(
			size_t queryIdx,
			EncodedFastaSequence & querySeq,
			size_t maxQueryLength,
			size_t subjectIdx,
			EncodedFastaSequence & subjectSeq,
			size_t maxSubjectLength,
			KnnHeap<Ranking> & knn,
			HausdorffCalculator * calculator
		) {
			double distance = calculator->ComputeDistance( *querySeq.base, *subjectSeq.base );

#if 0
#if USE_OMP
#pragma omp critical
#endif
			{
				cerr << querySeq->IdStr() << " " << subjectSeq->IdStr() << "\n";

				auto & rowMinima = calculator->GetRowMinima();
				cerr << "query.Sequence().size() = " << querySeq->Sequence().size() << "\n";
				cerr << "queryKmerCount          = " << calculator->queryKmerCount << "\n";
				cerr << "queryFragCount          = " << calculator->queryFragCount << "\n";
				cerr << "rowMinima               = ";

				double rowSum = 0;

				for ( uint i = 0; i < calculator->queryFragCount; i++ ) {
					cerr << "\t" << rowMinima[i];
					rowSum += rowMinima[i];
				}
				cerr << "\n";

				double rowAverage = rowSum / calculator->queryFragCount;
				cerr << "row average             = " << rowAverage << "\n";

				auto & colMinima = calculator->GetColMinima();
				cerr << "subject.Sequence().size() = " << subjectSeq->Sequence().size() << "\n";
				cerr << "subjectKmerCount          = " << calculator->subjectKmerCount << "\n";
				cerr << "subjectFragCount          = " << calculator->subjectFragCount << "\n";
				cerr << "colMinima                 = ";

				double colSum = 0;

				for ( uint i = 0; i < calculator->subjectFragCount; i++ ) {
					cerr << "\t" << colMinima[i];
					colSum += colMinima[i];
				}
				cerr << "\n";

				double colAverage = colSum / calculator->subjectFragCount;
				cerr << "col average             = " << colAverage << "\n";
				cerr << "distance                = " << (rowAverage + colAverage) / 2 << "\n";

				cerr << "\n\n-----------------------------\n";
			}
#endif

			Ranking r( querySeq.IdStr(), subjectSeq.IdStr(), distance, 0, 0 );
			knn.push( r );

#if WANT_DIAGNOSTIC_STREAM
			if ( diagnosticStream ) {
#if USE_OMP
#pragma omp critical
#endif
				{
					ostream & out( *diagnosticStream );

					out << querySeq.IdStr() << "," << subjectSeq.IdStr() << "," << distance << "," << 0 << "\n";

					out << querySeq.IdStr();

					for ( size_t i = 0; i < calculator->queryFragCount; i++ ) {
						int distance = calculator->rowMinima[i];
						out << "," << distance;
					}

					out << "\n" << subjectSeq.IdStr();

					for ( size_t i = 0; i < calculator->subjectFragCount; i++ ) {
						int distance = calculator->colMinima[i];
						out << "," << distance;
					}

					out << "\n" << querySeq.IdStr() << "," << subjectSeq.IdStr();

					for ( size_t i = 0; i < calculator->queryFragCount; i++ ) {
						out << ";";

						for ( size_t j = 0; j < calculator->subjectFragCount; j++ ) {
							if ( j > 0 ) {
								out << ",";
							}
							out << calculator->kmerDistCache( i, j );
						}
					}

					out << "\n";
				}
			}
#endif
		}
	};
}
