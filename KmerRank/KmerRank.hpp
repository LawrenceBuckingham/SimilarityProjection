#pragma once

#include "QutBio.h"
#include "Homologs.hpp"
#include "KmerDistanceCache.hpp"

#include <set>
#include <iomanip>
#include <numeric>
#include <cstdio>

using namespace QutBio;
using namespace std;

#undef max

// TODO: Use AAClustSigEncode to generate database and query sequences, and load signatures instead of the full cluster file to obtain posting lists.

namespace KmerQuery {
	class Main {

		using DistanceFunction = KmerDistanceCache2;
		using Cluster = KmerCluster<DistanceFunction, Kmer>;
		using Codebook = KmerCodebook<DistanceFunction, Kmer>;
		using Fasta = FastaSequence;
		using Seq = EncodedFastaSequence;
		using Proto = KmerClusterPrototype;

	public:
		/// <summary> Main class contains the program entry point: Run().
		/// </summary>
		static void Run( int argc, char ** argv ) {
			Args arguments( argc, argv );

			arguments.Show();

			if ( arguments.IsDefined( "numThreads" ) ) {
				int numThreads;
				if ( arguments.Get( "numThreads", numThreads ) ) {
					omp_set_num_threads( numThreads );
					cerr << "OMP thread count set to " << numThreads << endl;
				}
			}
			else {
				omp_set_num_threads( 1 );
			}

			try {
				auto start = omp_get_wtime();
				Process( arguments );
				auto end = omp_get_wtime();
				FILE *f = fopen( "KmerRank_time.txt", "a" );
				fprintf( f, "Elapsed time: %fs\n", (end - start) );
				fclose( f );
			}
			catch ( Exception & e ) {
				string message = e.what();

				if ( message.find( "Safe to ignore" ) == string::npos ) {
					(cerr << e.File() << "(" << e.Line() << "): " << message << "\n").flush();
					exit( 1 );
				}
			}
		}

	private:

		static void WriteRankingsToFile(
			vector<Ranking> & rankings,
			size_t maxRecords,
			ofstream & rankingFile
		) {
			const size_t N = std::min( rankings.size(), maxRecords );

			if ( N == 0 ) return;

			rankingFile << *(rankings[0].queryId) << "," << N;

			for ( size_t i = 0; i < N; i++ ) {
				// rankingFile << rankings[i] << '\n';
				rankingFile << "," << *(rankings[i].subjectId) << "," << (-rankings[i].distance);
			}

			rankingFile << "\n";
		}

		static void Process( Args & arguments ) {
			KmerSequenceRanker_Params parms{ arguments };

			cout << "--seed " << parms.seed << endl;

			using DistanceFunction = KmerDistanceCache2;
			using KmerType = Kmer;
			using Codebook = KmerCodebook<DistanceFunction, KmerType>;

			BlosumDifferenceFunction dist( parms.matrix );
			DistanceFunction distanceFunction( parms.alphabet, &dist );

			bool filesOk = true;

			if ( parms.skip > 1 && parms.fragLength > 1 ) {
				filesOk = false;
				(cerr << "Present version cannot manage skip > 1 and fragLength > 1 simultaneously.\n").flush();
			}

			if ( !File::Exists( parms.dbFile ) ) {
				filesOk = false;
				(cerr << "Database file " << parms.dbFile << " cannot be opened to read." << endl).flush();
			}

			if ( !File::Exists( parms.queryFile ) ) {
				filesOk = false;
				(cerr << "Query file " << parms.queryFile << " cannot be opened to read." << endl).flush();
			}

			if ( !filesOk ) {
				throw Exception( "Problem found with parameters. See stderr.", FileAndLine );
			}

			auto dbSeqs = Load::Fasta( parms.dbFile, parms.idIndex, parms.alphabet );
			PreprocessDataset( dbSeqs, parms.kmerLength, parms.alphabet->DefaultSymbol() );

			auto db = Load::Encoded( dbSeqs, parms.classIndex, parms.alphabet, parms.kmerLength, distanceFunction.CharsPerWord(), 'x' );
			PreprocessDataset( db );
			
			Index<Seq> seqIdx( db );
			KmerIndex kmerIdx( db, parms.kmerLength );

			(cerr << arguments.ProgName() << ": " << db.size() << " sequences loaded from '"
				<< parms.dbFile << "'.\n").flush();

			vector<FastaSequence *> querySeqs_;
			vector<Seq *> query_;

			if ( parms.queryFile != parms.dbFile ) {
				querySeqs_ = Load::Fasta( parms.queryFile, parms.idIndex, parms.alphabet );
				query_ = Load::Encoded( querySeqs_, parms.classIndex, parms.alphabet, parms.kmerLength, distanceFunction.CharsPerWord(), 'x' );
				PreprocessDataset( query_ );
				(cerr << arguments.ProgName() << ": Query dataset contains " << query_.size() << " sequences." << endl).flush();
			}

			vector<Seq *> &query = query_.size() > 0 ? query_ : db;
			vector<Seq *> querySubset;

			if ( parms.sampleSize <= 0 || parms.sampleSize >= query.size() ) {
				if ( parms.queryIdFile.size() > 0 ) {
					if ( !File::Exists( parms.queryIdFile ) ) {
						ostringstream cerr;
						cerr << arguments.ProgName() << ": " "queryIdFile '" << parms.queryIdFile << "' cannot be opened to read.\n";
						throw Exception( cerr.str(), FileAndLine );
					}

					ifstream queryIdStream( parms.queryIdFile );
					string queryId;

					while ( !(queryIdStream.eof() || queryIdStream.fail()) ) {
						getline( queryIdStream, queryId );

						if ( queryId.length() == 0 ) continue;

						try {
							querySubset.push_back( seqIdx.at( queryId ) );
						}
						catch ( std::out_of_range &ex ) {
							ostringstream cerr;
							cerr << arguments.ProgName() << ": " "Query sequence " << queryId << " not found in database.\n";
						}
					}
				}
				else {
					querySubset = query;
				}
			}
			else {
				UniformRealRandom rand( parms.seed );
				Selector want( rand, parms.sampleSize, query.size() );

				for ( auto seq : query ) {
					if ( want.SelectThis() ) {
						querySubset.push_back( seq );
					}
				}

				assert_equal( (size_t) parms.sampleSize, querySubset.size() );
			}

			// ( cerr << __FILE__ << ":" << __LINE__ << "\n" ).flush();

			(cerr << "Query subset contains " << querySubset.size() << " sequences." << endl).flush();

			// Write out the query IDs if requested.
			if ( parms.sampleSize > 0 && parms.sampleSize <= query.size() && parms.queryIdFile.length() > 0 ) {
				auto f = fopen( parms.queryIdFile.c_str(), "wb" );

				for ( auto qSeq_ : querySubset ) {
					auto & qSeq = *qSeq_;
					fprintf( f, "%s\n", qSeq.IdStr().c_str() );
				}

				fclose( f );
			}

			ofstream rankingFile( parms.rankingFile );

			rankingFile << "rankings," <<  querySubset.size() << "\n";

			Action1<vector<Ranking> &> postProcessRankings = [&]( vector<Ranking> & rankings ) {
				WriteRankingsToFile( rankings, parms.maxRecords, rankingFile );
			};

			Codebook * codebook = 0;
			vector<FastaSequence *> protoSeqs;
			vector<Proto *> prototypes;

			if ( parms.codebookFile.size() > 0 ) {
				if ( !File::Exists( parms.prototypeFile ) ) {
					cerr << "prototypeFile '" << parms.prototypeFile << "' does not exist.\n";
					throw Exception( "File not found.", FileAndLine );
				}

				protoSeqs = Load::Fasta( parms.prototypeFile, 0, parms.alphabet );
				prototypes = Load::Prototypes( protoSeqs, parms.alphabet, parms.kmerLength, distanceFunction.CharsPerWord() );
				Index<Proto> protoIndex( prototypes );

				FILE * f = fopen( parms.codebookFile.c_str(), "r" );

				if ( !f ) {
					ostringstream cerr;
					cerr << "Codebook " << parms.codebookFile << " cannot be opened for reading.\n";
					throw Exception( cerr.str(), FileAndLine );
				}

				codebook = new Codebook(
					parms.alphabet,
					distanceFunction,
					distanceFunction.CharsPerWord(),
					parms.kmerLength,
					seqIdx,
					protoIndex,
					kmerIdx,
					f,
					parms.balancedClusterSize
				);
				fclose( f );

				if ( codebook->Size() == 0 ) {
					ostringstream cerr;
					cerr << "Codebook contains no entries; run terminated.\n";
					throw Exception( cerr.str(), FileAndLine );
				}

				codebook->AllocateKmersToThreads( parms.numThreads );
			}
			else {
				cerr << "Information: argument --codebookFile not supplied. Exact Similarity Projection will be used.\n";
			}

#if defined(USE_BIG_CACHE)
			auto distance( RawKmerDistanceFunctionFactory::Factory( parms.distance, parms.matrixId ) );
			KmerDistanceCache cachedCalculator( parms.alphabet, distance );
#endif

			KmerSequenceRanker<DistanceFunction, KmerType> ranker(
				parms.matrix,
				parms.kmerMode,
				parms.fragMode,
				parms.alphabet,
				parms.fragLength,
				parms.kmerLength,
				parms.thresholdDistance,
				parms.defaultDistance,
				parms.pushKmerDistances,
				distanceFunction,
				parms.skip,
				parms.maxRecords
			);
			ranker.SetCodebook( codebook );
			ranker.SetQueryComplete( postProcessRankings );

#if defined(WANT_DIAGNOSTIC_STREAM)
			ostream * diagnosticStream = 0;

			if ( parms.diagnosticsFile.size() > 0 ) {
				cerr << "here!\n";
				diagnosticStream = new ofstream( parms.diagnosticsFile );
				ranker.SetDiagnosticStream( diagnosticStream );
			}
#endif
			ranker.RunJob( querySubset, db );
			rankingFile.close();

#if defined(WANT_DIAGNOSTIC_STREAM)
			if ( diagnosticStream ) {
				diagnosticStream->flush();
				delete diagnosticStream;
			}
#endif
		}

		/**
		*	<summary>
		*		If necessary, pads all sequences out to be as long as the word-size, k.
		*		Encodes the sequences to produce packed word matrices (if we insist that
		*		k be even, then we can collapse the matrices down to a single array, which
		*		would be nice.
		*	</summary>
		*/
		static void PreprocessDataset( vector<FastaSequence *> & db, size_t kmerLength, Symbol defaultSymbol ) {
#if USE_OMP
#pragma omp parallel for
#endif
			for ( size_t i = 0; i < db.size(); i++ ) {
				db[i]->EnsureLengthAtLeast( kmerLength, defaultSymbol );
			}
		}

		/**
		*	<summary>
		*		If necessary, pads all sequences out to be as long as the word-size, k.
		*		Encodes the sequences to produce packed word matrices (if we insist that
		*		k be even, then we can collapse the matrices down to a single array, which
		*		would be nice.
		*	</summary>
		*/
		static void PreprocessDataset( vector<Seq *> & db ) {
#if USE_OMP
#pragma omp parallel for
#endif
			for ( int i = 0; i < (int) db.size(); i++ ) {
				auto & seq = *db[i];
				seq.position = i;
			}
		}
	};

}


