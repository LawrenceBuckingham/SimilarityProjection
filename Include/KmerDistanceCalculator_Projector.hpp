#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif


#include "Alphabet.hpp"
#include "FastaSequence.hpp"
#include "FreeList.hpp"
#include "Fragment.hpp"
#include "KmerDistanceCache.hpp"
#include "KmerSequenceRanker_Params.hpp"
#include "Ranking.hpp"
#include "AveragePrecision.hpp"
#include "Util.hpp"
#include "Projector.hpp"
#include "BestHitDetector.hpp"

#include <string>
#include <ctime>
#include <cfloat>

#undef max
#undef min

namespace QutBio {

	class KmerDistanceCalculator_Projector: public KmerDistanceCalculator {
	protected:
		Alphabet * alphabet;
		SimilarityMatrix * matrix;

		uint kmerLength;

		const vector<EncodedFastaSequence *> * queryDataset;
		const vector<EncodedFastaSequence *> * subjectDataset;

		FragmentAggregationMode * kmerMode;

	public:
		/// <summary>
		/// Set this to a non-null stream to enable logging of kmer distances during calculations.
		/// </summary>
		ostream * fragDistOutStream = 0;

		/// <summary>
		///	Set to true to dump the kmer distance matrix to fragDistOutStream.
		/// </summary>
		bool dumpKmerDistCache = false;

		FreeList<Projector> projectors;

#if defined(USE_AMP_PROJECTOR)
		int charDistMatrix[1 << 16];
		array_view<const int, 2> charDistGPU;
#endif

		KmerDistanceCalculator(
			DistanceType * distanceType,
			uint matrixId,
			FragmentAggregationMode * kmerMode,
			FragmentAggregationMode * fragMode_ignored,
			Alphabet * alphabet,
			uint fragLength_ignored,
			uint interval_ignored,
			uint kmerLength
			) :
			matrix(GetMatrix(distanceType, matrixId)),
			alphabet(alphabet),
			kmerLength(kmerLength),
			kmerMode(kmerMode)
#if defined(USE_AMP_PROJECTOR)
			, charDistGPU {
			Projector::GetCharDistMatrix(matrix, charDistMatrix);
		}
#else
		{
		}
#endif

		virtual ~KmerDistanceCalculator() {}

		/// <summary>
		/// Called once before the start of all processing.
		/// </summary>

		virtual void setup(void) = 0;

		/// <summary>
		/// Called before each query is processed against the database.
		/// </summary>

		virtual void preProcess(void) = 0;

		/// <summary>
		/// Called to process each (query seq, db seq) pair.
		/// </summary>

		virtual void process(FastaSequence * querySeq, FastaSequence * subjectSeq, int subjectIdx, double distance) = 0;

		/// <summary>
		/// Called after all d b sequences have been processed for given query sequence.
		/// </summary>

		virtual void postProcess(void) {
			throw Exception("Call the other function!", FileAndLine);
		}

		virtual void postProcess(vector<int> & rowTotals, vector<int> & colTotals) {
			throw Exception("Call the other function!", FileAndLine);
		}

		/// <summary>
		/// Called once after all queries have been processed.
		/// </summary>

		virtual void cleanup(void) = 0;

		static SimilarityMatrix * GetMatrix(DistanceType * distanceType, uint matrixId) {
			return distanceType == DistanceType::HalperinEtAl() ? SimilarityMatrix::GetBlosum(matrixId)
				: distanceType == DistanceType::BlosumDistance() ? SimilarityMatrix::GetBlosum(matrixId)
				: 0;
		}

		/// Run the job.

		void RunJob(
			const vector<EncodedFastaSequence *> & query,
			const vector<EncodedFastaSequence *> & db,
			function<bool(const string &, const string &)> isHomolog
			) {
			this->queryDataset = &query;
			this->subjectDataset = &db;

			int n = (int) db.size();

			int maxQueryLength = numeric_limits<int>::min();

			for ( auto seq : query ) {
				int len = (int) seq->Sequence().size();

				if ( len > maxQueryLength ) maxQueryLength = len;
			}

			int maxSubjectLength = numeric_limits<int>::min();

			for ( auto seq : db ) {
				int len = (int) seq->Sequence().size();

				if ( len > maxSubjectLength ) maxSubjectLength = len;
			}

			int maxQueryKmerCount = maxQueryLength - kmerLength + 1;
			int maxSubjectKmerCount = maxSubjectLength - kmerLength + 1;

			setup( );

			for ( int queryIdx = 0; queryIdx < query.size(); queryIdx++ ) {
				auto querySeq = query[queryIdx];

#if defined(USE_AMP_PROJECTOR)
				vector<int> rowTotal(n);
				array_view<int, 1> rowTotalGPU(n, rowTotal);
				rowTotalGPU.discard_data();

				vector<int> colTotal(n);
				array_view<int, 1> colTotalGPU(n, colTotal);
				colTotalGPU.discard_data();
#endif

				preProcess();

				#pragma omp parallel for 
				for ( int subjectIdx = 0; subjectIdx < n; subjectIdx++ ) {
					EncodedFastaSequence * subjectSeq = db[subjectIdx];

#if defined(USE_AMP_PROJECTOR)
					auto projectorFactory = [=]() {
						return new Projector(matrix, kmerLength, charDistGPU, maxQueryKmerCount, maxSubjectKmerCount );
					};

					auto projector = projectors.Allocate( projectorFactory );
					projector->ComputeDistance(querySeq, subjectSeq, rowTotalGPU, colTotalGPU, subjectIdx);
					process(querySeq, subjectSeq, subjectIdx, 0);

					projectors.Free( projector );
#else
					SequenceDistanceFunction * f = kmerMode == FragmentAggregationMode::BestOfBest()
						? (SequenceDistanceFunction *) new BestHitDetector(matrix, kmerLength)
						: (SequenceDistanceFunction *) new Projector(matrix, kmerLength);

					double dist = f->ComputeDistance(querySeq, subjectSeq);
					process(querySeq, subjectSeq, subjectIdx, dist);

					delete f;
#endif
				}

#if defined(USE_AMP_PROJECTOR)
				rowTotalGPU.synchronize();
				colTotalGPU.synchronize();
				postProcess(rowTotal, colTotal);
#else
				postProcess();
#endif
			}

			cleanup();
		}

		const SimilarityMatrix & Matrix() {
			return *matrix;
		}

	private:
	};
#endif
	}
