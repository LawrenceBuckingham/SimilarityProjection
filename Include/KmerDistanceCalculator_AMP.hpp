#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "Alphabet.hpp"
#include "Console.hpp"
#include "Delegates.hpp"
#include "EncodedKmer.hpp"
#include "FastaSequence.hpp"
#include "Fragment.hpp"
#include "KmerDistanceCache.hpp"
#include "KmerSequenceRanker_Params.hpp"
#include "Ranking.hpp"
#include "AveragePrecision.hpp"
#include "Util.hpp"

#include <string>
#include <ctime>
#include <cfloat>

#include <amp.h>
using namespace concurrency;

#undef max
#undef min

namespace QutBio {
#if defined(USE_AMP)
	class KmerDistanceCalculator {
	protected:
		Alphabet * alphabet;
		SimilarityMatrix * matrix;

		// Copy of matrix->dict, formatted for use on GPU
		int charDist[1 << 16];

		FragmentAggregationMode * kmerMode;
		FragmentAggregationMode * fragMode;

		/// <summary>
		///		kmer distance cache with max(queryKmerCount) rows and max(subjectKmerCount) columns.
		///	</summary>
		int * kmerDistCache = 0;

		Array<Array<double>> fragDistCache;
		int vocabSize;

		uint fragmentLength;
		uint interval;
		uint kmerLength;

#if defined(USE_BIG_CACHE)
		size_t charsPerWord;
		KmerDistanceCache & cachedCalculator;
		KmerDistanceFunction kmerCodeDistance;
#endif // defined(USE_BIG_CACHE)


		using FragmentDistance = double( *)(
			KmerDistanceCalculator & instance,
			size_t queryStart,
			size_t subjectStart,
			size_t fragmentLength,
			double threshold
			);

		FragmentDistance fragmentDistance;

		using CollectionDistance = double( *)(
			KmerDistanceCalculator & instance,
			FragmentDistance fragDistance,
			double threshold
			);
		CollectionDistance collectionDistance;

		uint64_t totalFragDistCalcs = 0;
		uint64_t fragsEliminatedByDistance = 0;

		// Query sequence data.
		FastaSequence * querySeq;
		const char * queryChars;
		size_t queryIdx;
		size_t queryKmerCount;
		size_t queryFragCount;

		// Subject (DB) sequence data.
		FastaSequence * subjectSeq;
		const char * subjectChars;
		size_t subjectIdx;
		size_t subjectKmerCount;
		size_t subjectFragCount;

#if defined(USE_BIG_CACHE)
		Array<EncodedKmer> queryTiles;
		Array<EncodedKmer> subjectTiles;
#endif // defined(USE_BIG_CACHE)

		const vector<FastaSequence *> * queryDataset;
		const vector<FastaSequence *> * subjectDataset;

	public:
		/// <summary>
		/// Set this to a non-null stream to enable logging of kmer distances during calculations.
		/// </summary>
		ostream * fragDistOutStream = 0;

		/// <summary>
		///	Set to true to dump the kmer distance matrix to fragDistOutStream.
		/// </summary>
		bool dumpKmerDistCache = false;

		KmerDistanceCalculator(
			/**/	DistanceType * distanceType,
			/**/	uint matrixId,
			/**/	FragmentAggregationMode * kmerMode,
			/**/	FragmentAggregationMode * fragMode,
			/**/	Alphabet * alphabet,
			/**/	uint fragLength,
			/**/	uint interval,
#if defined(USE_BIG_CACHE)
			/**/	uint kmerLength,
			/**/	KmerDistanceCache & cachedCalculator
#else
			/**/	uint kmerLength
#endif // defined(USE_BIG_CACHE)
			) :
#if defined(USE_BIG_CACHE)
			/**/	cachedCalculator( cachedCalculator ),
#endif // defined(USE_BIG_CACHE)

			/**/
			/**/	matrix( distanceType == DistanceType::HalperinEtAl()
			/**/ ? SimilarityMatrix::GetBlosum( matrixId )
			/**/ : distanceType == DistanceType::BlosumDistance() ? SimilarityMatrix::GetBlosum( matrixId )
			/**/ : 0 ),
			/**/
			/**/	alphabet( alphabet ),
			/**/
			/**/	fragmentLength( fragLength ),
			/**/
			/**/	interval( interval ),
			/**/
			/**/	kmerLength( kmerLength ),
			/**/
			/**/	kmerMode( kmerMode ),
			/**/
			/**/	fragMode( fragMode )
			/**/
		{

#if defined(USE_BIG_CACHE)
			charsPerWord = cachedCalculator.CharsPerWord();
#endif // defined(USE_BIG_CACHE)


			FragmentDistance fragmentDistances[] = {
				KmerBestOfBest,
				KmerHausdorff,
				KmerHausdorffAverage,
				KmerHausdorffAverageAverage,
				KmerSlice,
				KmerSliceVertical,
				KmerSliceNoFollow,
				KmerSliceVerticalNoFollow,
			};

			fragmentDistance = fragmentDistances[kmerMode->Value()];

			CollectionDistance collectionDistances[] = {
				FragBestOfBest,
				FragHausdorff,
				FragHausdorffAverage,
				FragHausdorffAverageAverage
			};

			collectionDistance = collectionDistances[fragMode->Value()];

#if defined(USE_BIG_CACHE)
			if ( kmerLength <= 18 ) {
				auto kmerDistanceFuncs = cachedCalculator.KmerDistanceFunctions();
				kmerCodeDistance = kmerDistanceFuncs[kmerLength - 1];
			}
			else {
				kmerCodeDistance = cachedCalculator.GetKmerDistances_Any;
			}
#endif // defined(USE_BIG_CACHE)

			GetCharDistMatrix( matrix, charDist );
		}

		virtual ~KmerDistanceCalculator() {}

		/// <summary>
		/// Called once before the start of all processing.
		/// </summary>

		virtual void setup( void ) = 0;

		/// <summary>
		/// Called before each query is processed against the database.
		/// </summary>

		virtual void preProcess( void ) = 0;

		/// <summary>
		/// Called to process each (query seq, db seq) pair.
		/// </summary>

		virtual bool process( void ) = 0;

		/// <summary>
		/// Called after all d b sequences have been processed for given query sequence.
		/// </summary>

		virtual void postProcess( void ) = 0;

		/// <summary>
		/// Called once after all queries have been processed.
		/// </summary>

		virtual void cleanup( void ) = 0;

		/// Run the job.

		void RunJob(
			const vector<FastaSequence *> & query,
			const vector<FastaSequence *> & db,
			function<bool( const string &, const string & )> isHomolog
			) {
			PreallocateKmerDistCache( query, db );
			setup();

			bool keep_going = true;
			bool precomputeKmers = kmerMode->Value() < FragmentAggregationMode::Slice()->Value();

			for ( queryIdx = 0; keep_going && queryIdx < query.size(); queryIdx++ ) {
				querySeq = query[queryIdx];
				queryChars = querySeq->Sequence().c_str();
				queryKmerCount = querySeq->Sequence().size() - kmerLength + 1;
				queryFragCount = Fragment::GetCount( queryKmerCount, fragmentLength, interval );

#if defined(USE_BIG_CACHE)
				Fragment::DissectString( queryChars, querySeq->Sequence().size(), kmerLength, 1, 1,
					[&]( int fragmentIndex, int fragmentCount, const char * subSequence, size_t len ) {
					alphabet->Encode( subSequence, len, charsPerWord, queryTiles[fragmentIndex] );
				}
				);
#endif // defined(USE_BIG_CACHE)

				assert_true( queryFragCount <= fragDistCache.Length() );

				if ( fragDistOutStream ) {
					( *fragDistOutStream ) << "query," << querySeq->IdStr() << ",fragments," << queryFragCount << endl;
				}

				preProcess();

				for ( subjectIdx = 0; keep_going && subjectIdx < db.size(); subjectIdx++ ) {
					subjectSeq = db[subjectIdx];
					subjectChars = subjectSeq->Sequence().c_str();
					subjectKmerCount = subjectSeq->Sequence().size() - kmerLength + 1;
					subjectFragCount = Fragment::GetCount( subjectKmerCount, fragmentLength, interval );

#if defined(USE_BIG_CACHE)
					Fragment::DissectString( subjectChars, subjectSeq->Sequence().size(), kmerLength, 1, 1,
						[&]( int fragmentIndex, int fragmentCount, const char * subSequence, size_t len ) {
						alphabet->Encode( subSequence, len, charsPerWord, subjectTiles[fragmentIndex] );
					}
					);
#endif // defined(USE_BIG_CACHE)

					assert_true( subjectFragCount <= fragDistCache[0].Length() );

					if ( fragDistOutStream ) {
						( *fragDistOutStream ) << "subject," << db[subjectIdx]->IdStr() << ",fragments," << subjectFragCount << ",relevant," << ( isHomolog( querySeq->IdStr(), db[subjectIdx]->IdStr() ) ? "true" : "false" ) << endl;
					}

					if ( precomputeKmers ) {
						UpdateKmerDistCache();
					}

					keep_going = process();
				}

				postProcess();
			}

			cleanup();

			delete[] kmerDistCache;
		}

		const SimilarityMatrix & Matrix() {
			return *matrix;
		}


	private:
		size_t maxQueryLength;
		size_t maxSubjectLength;

		void PreallocateKmerDistCache(
			const vector<FastaSequence *> & query,
			const vector<FastaSequence *> & db
			) {
			queryDataset = &query;
			subjectDataset = &db;
			maxQueryLength = GetMaxLength( query );
			maxSubjectLength = GetMaxLength( db );

#if defined(USE_BIG_CACHE)
			queryTiles.Resize( maxQueryLength );
			subjectTiles.Resize( maxSubjectLength );
#endif // defined(USE_BIG_CACHE)

			kmerDistCache = new int[maxQueryLength * maxSubjectLength];

			int maxQueryFragCount = Fragment::GetCount( maxQueryLength, fragmentLength, interval );
			int maxSubjectFragCount = Fragment::GetCount( maxSubjectLength, fragmentLength, interval );

			fragDistCache.Resize( maxQueryFragCount );

			for ( int i = 0; i < maxQueryFragCount; i++ ) {
				fragDistCache[i].Resize( maxSubjectFragCount );
			}
		}

		static uint GetMaxLength( const vector<FastaSequence *> & dataset ) {
			uint queryMax = 0;

			for ( uint i = 0; i < dataset.size(); i++ ) {
				uint len = (uint) dataset[i]->Sequence().size();

				if ( len > queryMax ) {
					queryMax = len;
				}
			}

			return queryMax;
		}

		// Constructs a character mutual distance table from a similarity matrix. 
		static void GetCharDistMatrix(
			const SimilarityMatrix * m,
			int charDist[]
			) {
			for ( int i = 0; i < 128; i++ ) {
				for ( int j = 0; j < 128; j++ ) {
					int distance = m->Difference( (char) i, (char) j );
					charDist[i * 128 + j] = distance;
				}
			}
		}

		/// <summary>
		///	Given strings s and t, having length m and n respectively, and a similarity matrix 
		///	reformatted as a 128 by 128 array of integers, populates a vector contain a kmer mutual
		///	distance table.
		/// </summary>
		///	<param name="queryBytes"></param>
		static void GetKmerDistanceMatrix(
			const char * queryChars,
			const char * subjectChars,
			int k,
			const int charDist[128 * 128],
			int * distance
			) {
			const int m = (int) strlen( queryChars ) - k + 1;
			const int n = (int) strlen( subjectChars ) - k + 1;

			vector<int> sv( m + k - 1 ); for ( int i = 0; i < m + k - 1; i++ ) { sv[i] = queryChars[i]; }
			vector<int> tv( n + k - 1 ); for ( int i = 0; i < n + k - 1; i++ ) { tv[i] = subjectChars[i]; }

			array_view<const int, 1> a( m + k - 1, sv );
			array_view<const int, 1> b( n + k - 1, tv );
			array_view<const int, 2> cd( 128, 128, charDist );
			array_view<int, 2> c( m, n, distance );
			c.discard_data();

			parallel_for_each( c.extent,
				[=]( index<2> idx ) restrict( amp ) {
				int row = idx[0];
				int col = idx[1];
				int sum = 0;
				for ( int i = 0; i < k; i++ ) {
					int x = a[row + i];
					int y = b[col + i];
					sum += cd[x][y];
				}
				c[idx] = sum;
			}
			);

			c.synchronize();
		}

		// Compute the kmer mutual distance matrix.
		void UpdateKmerDistCache() {
			GetKmerDistanceMatrix( queryChars, subjectChars, kmerLength, charDist, kmerDistCache );
		}

		///	<summary>
		///		Calculates the threshold Hausdorff average (bidirectional) 
		///		distance between two fragments using the maximum of the two
		///		one way distances.
		///	</summary>
		static double FragHausdorffAverage(
			KmerDistanceCalculator & instance,
			FragmentDistance getFragmentDistance,
			double threshold
			) {
			double totXY = 0;
			int obsXY = 0;

			double totYX = 0;
			int obsYX = 0;

			double qRatio = ( (double) instance.queryKmerCount - instance.fragmentLength ) / ( instance.queryFragCount - 1 );
			double sRatio = ( (double) instance.subjectKmerCount - instance.fragmentLength ) / ( instance.subjectFragCount - 1 );

			// Compute one-way distance from x to y, and also cache distances.
			for ( uint i = 0; i < instance.queryFragCount; i++ ) {
				uint qStart = (uint) round( i * qRatio );
				double min = DBL_MAX;
				uint minIdx = numeric_limits<uint>::max();

				for ( uint j = 0; j < instance.subjectFragCount; j++ ) {
					uint sStart = (uint) round( j * sRatio );
					double distance = getFragmentDistance( instance, qStart, sStart, instance.fragmentLength, threshold );

					instance.fragDistCache[i][j] = distance;

					instance.totalFragDistCalcs++;

					if ( distance >= threshold ) {
						instance.fragsEliminatedByDistance++;
					}

					if ( distance < min ) {
						min = distance;
					}
				}

				if ( instance.fragDistOutStream ) {
					( *instance.fragDistOutStream )
						<< instance.querySeq->IdStr()
						<< ","
						<< i
						<< ","
						<< instance.subjectSeq->IdStr()
						<< ","
						<< minIdx
						<< ","
						<< min
						<< endl;
				}

				if ( min < threshold ) {
					obsXY++;
					totXY += min;
				}
			}

			// Compute one-way distance from y to x, using cache distances.
			for ( uint j = 0; j < instance.subjectFragCount; j++ ) {
				double min = DBL_MAX;
				uint minIdx = numeric_limits<uint>::max();

				for ( uint i = 0; i < instance.queryFragCount; i++ ) {
					double distance = instance.fragDistCache[i][j];

					if ( distance < min ) {
						min = distance;
						minIdx = i;
					}
				}

				if ( instance.fragDistOutStream ) {
					( *instance.fragDistOutStream )
						<< instance.subjectSeq->IdStr()
						<< ","
						<< j
						<< ","
						<< instance.querySeq->IdStr()
						<< ","
						<< minIdx
						<< ","
						<< min
						<< endl;
				}

				if ( min < threshold ) {
					obsYX++;
					totYX += min;
				}
			}

			double avgXY = obsXY > 0 ? totXY / obsXY : threshold;
			double avgYX = obsYX > 0 ? totYX / obsYX : threshold;

			double result = avgXY > avgYX ? avgXY : avgYX;

			if ( instance.fragDistOutStream ) {
				( *instance.fragDistOutStream )
					<< instance.querySeq->IdStr()
					<< ","
					<< instance.subjectSeq->IdStr()
					<< ","
					<< avgXY
					<< ","
					<< avgYX
					<< ","
					<< result
					<< endl;
			}

			return result;
		}

		///	<summary>
		///		Calculates the threshold Hausdorff average (bidirectional) 
		///		distance between two fragments using the average of the two 
		///		one-way fragment distances.
		///	</summary>
		static double FragHausdorffAverageAverage(
			KmerDistanceCalculator & instance,
			FragmentDistance getFragmentDistance,
			double threshold
			) {
			double totXY = 0;
			int obsXY = 0;

			double totYX = 0;
			int obsYX = 0;

			double qRatio = ( (double) instance.queryKmerCount - instance.fragmentLength ) / ( instance.queryFragCount - 1 );
			double sRatio = ( (double) instance.subjectKmerCount - instance.fragmentLength ) / ( instance.subjectFragCount - 1 );

			// Compute one-way distance from x to y, and also cache distances.
			for ( uint i = 0; i < instance.queryFragCount; i++ ) {
				uint qStart = (uint) round( i * qRatio );
				double min = DBL_MAX;
				uint minIdx = numeric_limits<uint>::max();

#pragma omp parallel for
				for ( int j = 0; j < (int) instance.subjectFragCount; j++ ) {
					uint sStart = (uint) round( j * sRatio );
					double distance = getFragmentDistance( instance, qStart, sStart, instance.fragmentLength, threshold );
					instance.fragDistCache[i][j] = distance;
					instance.totalFragDistCalcs++;
					if ( distance >= threshold ) {
						instance.fragsEliminatedByDistance++;
					}
				}

				for ( int j = 0; j < (int) instance.subjectFragCount; j++ ) {
					if ( instance.fragDistCache[i][j] < min ) {
						min = instance.fragDistCache[i][j];
						minIdx = j;
					}
				}

				if ( instance.fragDistOutStream ) {
					( *instance.fragDistOutStream )
						<< instance.querySeq->IdStr()
						<< ","
						<< i
						<< ","
						<< instance.subjectSeq->IdStr()
						<< ","
						<< minIdx
						<< ","
						<< min
						<< endl;
				}

				if ( min < threshold ) {
					obsXY++;
					totXY += min;
				}
			}

			// Compute one-way distance from y to x, using cache distances.
			for ( uint j = 0; j < instance.subjectFragCount; j++ ) {
				double min = DBL_MAX;
				uint minIdx = numeric_limits<uint>::max();

				for ( int i = 0; i < (int) instance.queryFragCount; i++ ) {
					double distance = instance.fragDistCache[i][j];
					if ( distance < min ) {
						min = distance;
						minIdx = i;
					}
				}

				if ( instance.fragDistOutStream ) {
					( *instance.fragDistOutStream )
						<< instance.subjectSeq->IdStr()
						<< ","
						<< j
						<< ","
						<< instance.querySeq->IdStr()
						<< ","
						<< minIdx
						<< ","
						<< min
						<< endl;
				}

				if ( min < threshold ) {
					obsYX++;
					totYX += min;
				}
			}

			double avgXY = obsXY > 0 ? totXY / obsXY : 1000000;
			double avgYX = obsYX > 0 ? totYX / obsYX : 1000000;

			double result = ( avgXY + avgYX ) / 2;

			if ( instance.fragDistOutStream ) {
				ostream & f = ( *instance.fragDistOutStream );

				f << instance.querySeq->IdStr()
					<< ","
					<< instance.subjectSeq->IdStr()
					<< ","
					<< avgXY
					<< ","
					<< avgYX
					<< ","
					<< result
					<< endl
					<< endl
					<< "Fragment distance Matrix:"
					<< endl;

				for ( uint j = 0; j < instance.subjectFragCount; j++ ) {
					f << "," << j;
				}

				f << endl;

				for ( uint i = 0; i < instance.queryFragCount; i++ ) {
					f << i;

					for ( uint j = 0; j < instance.subjectFragCount; j++ ) {
						f << "," << instance.fragDistCache[i][j];
					}

					f << endl;
				}

				f << endl;
			}

			return result;
		}

		static double FragHausdorff( KmerDistanceCalculator & instance, FragmentDistance getFragmentDistance, double threshold ) {
			double maxXY = -DBL_MAX;
			double maxYX = -DBL_MAX;

			double qRatio = ( (double) instance.queryKmerCount - instance.fragmentLength ) / ( instance.queryFragCount - 1 );
			double sRatio = ( (double) instance.subjectKmerCount - instance.fragmentLength ) / ( instance.subjectFragCount - 1 );

#pragma omp parallel for
			// Compute fragment distances and cache them.
			for ( int i = 0; i < (int) instance.queryFragCount; i++ ) {
				uint qStart = (uint) round( i * qRatio );

				for ( int j = 0; j < (int) instance.subjectFragCount; j++ ) {
					uint sStart = (uint) round( j * sRatio );
					double distance = getFragmentDistance( instance, qStart, sStart, instance.fragmentLength - instance.kmerLength + 1, threshold );
					instance.fragDistCache[i][j] = distance;
				}
			}

			// Compute one-way distance from x to y, and also cache distances.
			for ( uint qStart = 0, i = 0; i < instance.queryFragCount; qStart += instance.interval, i++ ) {
				double min = DBL_MAX;

				for ( uint sStart = 0, j = 0; j < instance.subjectFragCount; sStart += instance.interval, j++ ) {
					double distance = instance.fragDistCache[i][j];

					if ( distance < min ) {
						min = distance;
					}
				}

				if ( min > maxXY ) {
					maxXY = min;
				}
			}

			// Compute one-way distance from y to x, using cache distances.
			for ( uint sStart = 0, j = 0; j < instance.subjectFragCount; sStart += instance.interval, j++ ) {
				double min = DBL_MAX;

				for ( uint qStart = 0, i = 0; i < instance.queryFragCount; qStart += instance.interval, i++ ) {
					double distance = instance.fragDistCache[i][j];

					if ( distance < min ) {
						min = distance;
					}
				}

				if ( min > maxYX ) {
					maxYX = min;
				}
			}

			return maxXY > maxYX ? maxXY : maxYX;
		}

		static double FragBestOfBest(
			KmerDistanceCalculator &instance,
			FragmentDistance getFragmentDistance,
			double threshold
			) {
			double minDist = numeric_limits<double>::max();

			double qRatio = ( (double) instance.queryKmerCount - instance.fragmentLength ) / ( instance.queryFragCount - 1 );
			double sRatio = ( (double) instance.subjectKmerCount - instance.fragmentLength ) / ( instance.subjectFragCount - 1 );

#pragma omp parallel for
			// Compute fragment distances and cache them.
			for ( int i = 0; i < (int) instance.queryFragCount; i++ ) {
				uint qStart = (uint) round( i * qRatio );

				for ( int j = 0; j < (int) instance.subjectFragCount; j++ ) {
					uint sStart = (uint) round( j * sRatio );
					double distance = getFragmentDistance( instance, qStart, sStart, instance.fragmentLength - instance.kmerLength + 1, threshold );
					instance.fragDistCache[i][j] = distance;
				}
			}

			for ( uint qStart = 0, i = 0; i < instance.queryFragCount; qStart += instance.interval, i++ ) {
				for ( uint sStart = 0, j = 0; j < instance.subjectFragCount; sStart += instance.interval, j++ ) {
					double distance = instance.fragDistCache[i][j];

					if ( distance < minDist ) {
						minDist = distance;
					}
				}
			}

			return minDist;
		}

		/// <summary> Returns the min-of-min kmer distance between two sets of encoded kmers,
		/// </summary>
		/// <param name="qStart">The index of the first kmer in the query fragment.</param>
		/// <param name="sStart">The index of the first kmer in the reference fragment.</param>
		/// <param name="len">The number of kmers in the fragment.</param>
		///	<param name="threshold">Early stopping criterion for summative, non-negative, distances. 
		///		If cumulative distance equals or exceeds threshold, stop calculating and return cumulative 
		///		value at the stopping time.
		///	</param>
		/// <returns></returns>

		static double KmerBestOfBest(
			KmerDistanceCalculator & instance,
			size_t qStart,
			size_t sStart,
			size_t len,
			double threshold
			) {
			int minDist = MAX_DIST;
			size_t qEnd = min( qStart + len, instance.queryKmerCount );
			size_t sEnd = min( sStart + len, instance.subjectKmerCount );

			for ( size_t i = qStart; i < qEnd; i++ ) {
				int * row = instance.kmerDistCache + i*instance.subjectKmerCount;

				for ( size_t j = sStart; j < sEnd; j++ ) {
					int distance = row[j];

					if ( distance < minDist ) {
						minDist = distance;
					}
				}
			}

			return minDist;
		}

		static double KmerSlice(
			KmerDistanceCalculator & instance,
			size_t qStart,
			size_t sStart,
			size_t len,
			double threshold
			) {
			Distance minDist = MAX_DIST;
			size_t minI = 0;
			size_t minJ = 0;
			size_t qEnd = min( qStart + len, instance.queryKmerCount );
			size_t sEnd = min( sStart + len, instance.subjectKmerCount );

			for ( size_t i = qStart; i < qEnd; i++ ) {
				size_t j0 = ( sEnd - 1 ) - ( i - qStart ) * ( sEnd - sStart ) / ( qEnd - qStart );

				for ( int j = (int) j0; j >= (int) sStart && j >= (int) j0 - 1; j-- ) {
					// assert_true((int) sStart <= j && j < (int) sEnd);
					Distance distance = instance.matrix->Difference( instance.queryChars + i, instance.subjectChars + j, instance.kmerLength );

					if ( distance < minDist ) {
						minDist = distance;
						minI = i;
						minJ = j;
					}
				}
			}

			for ( size_t i = minI + 1, j = minJ + 1; i < qEnd && j < sEnd; i++, j++ ) {
				Distance distance = instance.matrix->Difference( instance.queryChars + i, instance.subjectChars + j, instance.kmerLength );

				if ( distance < minDist ) {
					minDist = distance;
				}
			}

			for ( int i = (int) minI - 1, j = (int) minJ - 1; i >= (int) qStart && j >= (int) sStart; i--, j-- ) {
				Distance distance = instance.matrix->Difference( instance.queryChars + i, instance.subjectChars + j, instance.kmerLength );

				if ( distance < minDist ) {
					minDist = distance;
				}
			}

			return minDist;
		}

		static double KmerSliceVertical(
			KmerDistanceCalculator & instance,
			size_t qStart,
			size_t sStart,
			size_t len,
			double threshold
			) {
			Distance minDist = MAX_DIST;
			size_t minI = 0;
			size_t minJ = 0;
			size_t qEnd = min( qStart + len, instance.queryKmerCount );
			size_t sEnd = min( sStart + len, instance.subjectKmerCount );

			for ( size_t i = qStart; i < qEnd; i++ ) {
				size_t j = sStart + ( sEnd - sStart ) / 2;

				Distance distance = instance.matrix->Difference(
					instance.queryChars + i, instance.subjectChars + j, instance.kmerLength
					);

				if ( distance < minDist ) {
					minDist = distance;
					minI = i;
					minJ = j;
				}
			}

			for ( size_t i = minI + 1, j = minJ + 1; i < qEnd && j < sEnd; i++, j++ ) {
				Distance distance = instance.matrix->Difference(
					instance.queryChars + i, instance.subjectChars + j, instance.kmerLength
					);

				if ( distance < minDist ) {
					minDist = distance;
				}
			}

			for ( int i = (int) minI - 1, j = (int) minJ - 1; i >= (int) qStart && j >= (int) sStart; i--, j-- ) {
				Distance distance = instance.matrix->Difference(
					instance.queryChars + i, instance.subjectChars + j, instance.kmerLength
					);

				if ( distance < minDist ) {
					minDist = distance;
				}
			}

			return minDist;
		}

		static double KmerSliceNoFollow(
			KmerDistanceCalculator & instance,
			size_t qStart,
			size_t sStart,
			size_t len,
			double threshold
			) {
			Distance minDist = MAX_DIST;
			size_t minI = 0;
			size_t minJ = 0;
			size_t qEnd = min( qStart + len, instance.queryKmerCount );
			size_t sEnd = min( sStart + len, instance.subjectKmerCount );

			for ( size_t i = qStart; i < qEnd; i++ ) {
				size_t j0 = ( sEnd - 1 ) - ( i - qStart ) * ( sEnd - sStart ) / ( qEnd - qStart );

				for ( int j = (int) j0; j >= (int) sStart && j >= (int) j0 - 1; j-- ) {
					assert_true( (int) sStart <= j && j < (int) sEnd );

					Distance distance = instance.matrix->Difference(
						instance.queryChars + i, instance.subjectChars + j, instance.kmerLength
						);

					if ( distance < minDist ) {
						minDist = distance;
						minI = i;
						minJ = j;
					}
				}
			}

			return minDist;
		}

		static double KmerSliceVerticalNoFollow(
			KmerDistanceCalculator & instance,
			size_t qStart,
			size_t sStart,
			size_t len,
			double threshold
			) {
			Distance minDist = MAX_DIST;
			size_t minI = 0;
			size_t minJ = 0;
			size_t qEnd = min( qStart + len, instance.queryKmerCount );
			size_t sEnd = min( sStart + len, instance.subjectKmerCount );

			for ( size_t i = qStart; i < qEnd; i++ ) {
				size_t j = sStart + ( sEnd - sStart ) / 2;

				Distance distance = instance.matrix->Difference(
					instance.queryChars + i, instance.subjectChars + j, instance.kmerLength
					);

				if ( distance < minDist ) {
					minDist = distance;
					minI = i;
					minJ = j;
				}
			}

			return minDist;
		}

		/// <summary> Returns the max(max-of-min) kmer distance between two sets of encoded kmers,
		/// </summary>
		/// <param name="qStart">The index of the first kmer in the query fragment.</param>
		/// <param name="sStart">The index of the first kmer in the reference fragment.</param>
		/// <param name="len">The number of kmers in the fragment.</param>
		///	<param name="threshold">Early stopping criterion for summative, non-negative, distances. 
		///		If cumulative distance equals or exceeds threshold, stop calculating and return cumulative 
		///		value at the stopping time.
		///	</param>
		/// <returns></returns>

		static double KmerHausdorff(
			KmerDistanceCalculator & instance,
			size_t qStart,
			size_t sStart,
			size_t len,
			double threshold
			) {
			int maxXY = numeric_limits<short>::min();
			int maxYX = numeric_limits<short>::min();
			size_t qEnd = min( qStart + len, instance.queryKmerCount );
			size_t sEnd = min( sStart + len, instance.subjectKmerCount );

			// Compute one-way distance from x to y, and also cache distances.
			for ( size_t i = qStart; i < qEnd; i++ ) {
				Distance min = numeric_limits<Distance>::max();
				int * row = instance.kmerDistCache + i * instance.subjectKmerCount;

				for ( size_t j = sStart; j < sEnd; j++ ) {
					Distance distance = row[j];

					if ( distance < min ) {
						min = distance;
					}
				}

				if ( min > maxXY ) {
					maxXY = min;
				}
			}

			// Compute one-way distance from y to x, using cache distances.
			for ( size_t j = sStart; j < sEnd; j++ ) {
				int min = numeric_limits<Distance>::max();

				for ( size_t i = qStart; i < qEnd; i++ ) {
					int distance = instance.kmerDistCache[i * instance.subjectKmerCount + j];

					if ( distance < min ) {
						min = distance;
					}
				}

				if ( min > maxYX ) {
					maxYX = min;
				}
			}

			return maxXY > maxYX ? maxXY : maxYX;
		}

		/// <summary> Returns the max(average-of-min) kmer distance between two sets of encoded kmers,
		/// </summary>
		/// <param name="qStart">The index of the first kmer in the query fragment.</param>
		/// <param name="sStart">The index of the first kmer in the reference fragment.</param>
		/// <param name="len">The number of kmers in the fragment.</param>
		/// <returns></returns>

		static double KmerHausdorffAverage(
			KmerDistanceCalculator & instance,
			size_t qStart,
			size_t sStart,
			size_t len,
			double threshold
			) {
			int totXY = 0;
			int totYX = 0;
			size_t qEnd = min( qStart + len, instance.queryKmerCount );
			size_t sEnd = min( sStart + len, instance.subjectKmerCount );

			int cutoffXY = (int) ( threshold * ( qEnd - qStart ) );
			int cutoffYX = (int) ( threshold * ( sEnd - sStart ) );

			// Compute one-way distance from x to y, using pre-cached cache distances.
			for ( size_t i = qStart; i < qEnd; i++ ) {
				int min = numeric_limits<Distance>::max();
				int * row = instance.kmerDistCache + i * instance.subjectKmerCount;

				for ( size_t j = sStart; j < sEnd; j++ ) {
					Distance distance = row[j];

					if ( distance < min ) {
						min = distance;
					}
				}

				totXY += min;

				// Early stop if average is going to exceed the threshold
				if ( totXY > cutoffXY ) {
					return (double) totXY / ( qEnd - qStart );
				}
			}

			// Compute one-way distance from y to x, using cache distances.
			for ( size_t j = sStart; j < sEnd; j++ ) {
				int min = numeric_limits<Distance>::max();

				for ( size_t i = qStart; i < qEnd; i++ ) {
					int distance = instance.kmerDistCache[i * instance.subjectKmerCount + j];

					if ( distance < min ) {
						min = distance;
					}
				}

				totYX += min;

				// Early stop if average is going to exceed the threshold
				if ( totYX > cutoffYX ) {
					return (double) totYX / ( sEnd - sStart );
				}
			}

			double avgXY = (double) totXY / ( qEnd - qStart );
			double avgYX = (double) totYX / ( sEnd - sStart );

			return avgXY > avgYX ? avgXY : avgYX;
		}

		/// <summary> Returns the average(average-of-min) kmer distance between two sets of encoded kmers,
		/// </summary>
		/// <param name="qStart">The index of the first kmer in the query fragment.</param>
		/// <param name="sStart">The index of the first kmer in the reference fragment.</param>
		/// <param name="len">The number of kmers in the fragment.</param>
		/// <returns></returns>

		static double KmerHausdorffAverageAverage(
			KmerDistanceCalculator & instance,
			size_t qStart,
			size_t sStart,
			size_t len,
			double threshold
			) {
			int totXY = 0;
			int totYX = 0;
			size_t qEnd = min( qStart + len, instance.queryKmerCount );
			size_t sEnd = min( sStart + len, instance.subjectKmerCount );

			// Compute one-way distance from x to y, using pre-cached cache distances.
			for ( size_t i = qStart; i < qEnd; i++ ) {
				int min = numeric_limits<Distance>::max();
				int * row = instance.kmerDistCache + i * instance.subjectKmerCount;

				for ( size_t j = sStart; j < sEnd; j++ ) {
					Distance distance = row[j];

					if ( distance < min ) {
						min = distance;
					}
				}

				totXY += min;
			}

			// Compute one-way distance from y to x, using cache distances.
			for ( size_t j = sStart; j < sEnd; j++ ) {
				Distance min = numeric_limits<Distance>::max();

				for ( size_t i = qStart; i < qEnd; i++ ) {
					int distance = instance.kmerDistCache[i * instance.subjectKmerCount + j];

					if ( distance < min ) {
						min = distance;
					}
				}

				totYX += min;
			}

			double avgXY = (double) totXY / ( qEnd - qStart );
			double avgYX = (double) totYX / ( sEnd - sStart );

			return ( avgXY + avgYX ) / 2;
		}

	};
#endif
}
