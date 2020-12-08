#include "Args.hpp"
#include "Assert.hpp"
#include "Delegates.hpp"
#include "Random.hpp"
#include "OmpTimer.h"
#include "kNearestNeighbours.hpp"
#include "Ranking.hpp"
#include "FragmentAggregationMode.hpp"
#include "Simproj.hpp"
#include "Sequence.hpp"
#include "Fragment.hpp"

#include <cstdio>
#include <stdlib.h>
#include <omp.h>
#include <set>
#include <vector>

using namespace QutBio;
using namespace std;

// Singletons.
Args* arguments;
mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

struct AASP {
public:

	struct Params {
	public:
		string dbFile;
		string queryFile;
		string outFile;
		size_t numThreads = 7;
		size_t wordLength = 0;
		int idIndex = 0;
		bool ok = true;
		int maxResults = 500;
		size_t fragLength = 1;
		FragmentAggregationMode* fragMode = FragmentAggregationMode::HausdorffAverageAverage();
		const SubstitutionMatrix* matrix = &SubstitutionMatrix::Blosum62();

		Params() {

			if ( arguments->IsDefined( "help" ) ) {
				vector<string> text{
"AASP:          Uses basic Similarity Projection with optional fragmentation to ",
"               get the top maxResults matches for a list of queries from a ",
"               database.",
"--help         Gets this text.",
"--dbFile       Required. A file path. The file will be parsed as a FASTA file ",
"               which contains amino acid sequences that have been clustered.",
"--queryFile    Required. The name of a file containing the prototypes."
"--outFile      Required. The name of a file which will be overwritten with ",
"               ranking records.",
"--idIndex      Required. The 0-origin position of the sequence ID field in ",
"               the pipe-separated definition line.",
"--wordLength   Required; The word length used for kmer tiling.",
"--numThreads   Optional; default value = 7. The number of OpenMP threads to ",
"               use in parallel regions.",
"--maxResults   Optional; default value = 500. The maximum number of rankings ",
"               to emit per query.",
"--fragMode     Optional. Fragment aggregation mode: one of the following:",
"               BestOfBest - the overall minimum fragment distance.",
"               Hausdorff - the maximum of the minimum one-way fragment ",
"               distances.",
"               HausdorffAverage - the maximum of the average one-way minimum ",
"               fragment distances.",
"               HausedorffAverageAverage - the average of the average one-way ",
"               minimum fragment distances.",
"               Default value is HausdorffAverageAverage.",
"--fragLength   Optional length of fragment. Default value = 1, yielding maximum ",
"				resolution. Within a fragment pair, the minimal pairwise word ",
"				distance is taken to be the distance between fragments. If frafLength",
"				is greater than 1, some fragments will be stretched or squashed to ",
"				ensure that the entire sequence is covered as evenly as possible.",
"--matrixFile   Optional path to text document containing a substitution matrix. ",
"				Default is built-in BLOSUM-62."
"               ",
				};

				for ( auto s : text ) {
					if ( s[0] == '-' && s[1] == '-' ) {
						cout << "\n";
					}
					cout << s << "\n";
				}
			}

			if ( !arguments->Get( "dbFile", this->dbFile ) ) {
				cout << arguments->ProgName() << ": error - required argument '--dbFile' not supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "queryFile", this->queryFile ) ) {
				cout << arguments->ProgName() << ": error - required argument '--queryFile' not "
					"supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "idIndex", idIndex ) ) {
				cout << arguments->ProgName() << ": error - required argument '--idIndex' not "
					"supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "numThreads", numThreads ) ) {
				cout << arguments->ProgName() << ": note - optional argument '--numThreads' not set"
					". Default value " << numThreads << " will be used.\n";
			}

			if ( !arguments->Get( "maxResults", maxResults ) ) {
				cout << arguments->ProgName() << ": note - optional argument '--maxResults' not set"
					". Default value " << maxResults << " will be used.\n";
			}

			if ( !arguments->Get( "wordLength", wordLength ) ) {
				cout << arguments->ProgName() << ": error: required argument '--wordLength' not supplied.\n";
				ok = false;
			}

			if ( !arguments->Get( "outFile", outFile ) ) {
				cout << arguments->ProgName() << ": note - required argument '--outFile' not provided.\n";
				ok = false;
			}

			if ( !arguments->Get( "fragMode", FragmentAggregationMode::Values(), fragMode ) ) {
				cout << arguments->ProgName() << ": note - optional argument '--fragMode' not provided.\n"
					<< "Default value " << fragMode->Name() << " will be used.\n";
			}

			if ( !arguments->Get( "fragLength", fragLength ) ) {
				cout << arguments->ProgName() << ": note - optional argument '--fragLength' not provided.\n"
					<< "Default value " << fragLength << " will be used.\n";
			}

			string matrixFile;
			if ( arguments->Get( "matrixFile", matrixFile ) ) {
				ifstream f(matrixFile);
				matrix = new SubstitutionMatrix(f);
			}
			else {
				cout << arguments->ProgName() << ": note - optional argument '--matrixFile' not provided.\n"
					<< "Built-in BLOSUM-62 will be used.\n";
			}
		}
	};

	static int Run() {
		Params p;

		if ( !p.ok ) {
			cout << "For help: " << arguments->ProgName() << " --help.\n";
			return 1;
		}

		omp_set_num_threads( p.numThreads );

		vector<Sequence*> db = Load( *p.matrix, p.dbFile, p.idIndex );;
		cout << arguments->ProgName() << ": " << db.size() << " reference sequences loaded.\n";

		vector<Sequence*> query = Load(  *p.matrix, p.queryFile, p.idIndex );
		cout << arguments->ProgName() << ": " << query.size() << " query sequences loaded.\n";

		auto cleanup = [&db, &query]() {
			Util::Free( query );
			Util::Free( db );
		};

		OMP_TIMER_DECLARE( rankTime );
		OMP_TIMER_START( rankTime );

		if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
			Rank<Simproj::BestOfBest>( query, db, p.wordLength, p.fragLength, *p.matrix, p.maxResults, p.outFile );

		else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
			Rank<Simproj::HausdorffAverageAverage>( query, db, p.wordLength, p.fragLength, *p.matrix, p.maxResults, p.outFile );

		else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
			Rank<Simproj::HausdorffAverage>( query, db, p.wordLength, p.fragLength, *p.matrix, p.maxResults, p.outFile );

		else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
			Rank<Simproj::Hausdorff>( query, db, p.wordLength, p.fragLength, *p.matrix, p.maxResults, p.outFile );

		else {
			cout << "Unknown fragMode.\n";
			cleanup();
			return 1;
		}


		OMP_TIMER_END( rankTime );

#if USE_OMP
		cout << "Ranking completed in " << OMP_TIMER( rankTime ) << "s.\n";
#endif

		cleanup();
		return 0;
	}

	static vector<Sequence *> Load(const SubstitutionMatrix & matrix, const string & fileName, int idIndex ) {
		ifstream f(fileName.c_str());
		LineReader r(f);
		vector<Sequence *> db = Sequence::ParseAllFasta( r, matrix, idIndex );
		return db;
	}

	/*
			size_t kmerCount = seq.KmerCount( parms.wordLength );
			size_t fragCount = Fragment::GetCount( kmerCount, parms.fragLength );
			double stepSize = Fragment::GetRealStepSize( kmerCount, parms.fragLength, fragCount );

			frags.resize( fragCount );

			auto process = [&frags, &summaryTfv, stepSize, alphabet]( EncodedFastaSequence* seq, size_t pos, size_t length ) {
				size_t fragIdx = (size_t) ( pos / stepSize );
				TermFreqVector& tfv = frags[fragIdx];
				Substring s( seq->Sequence().data(), pos, length, alphabet );
				tfv[s] ++;
				summaryTfv[s] ++;
			};

	*/

	template<double ( *Agg )( const vector<int>&, size_t, const vector<int>&, size_t )>
	static void Rank(
		const vector<Sequence*>& query,
		const vector<Sequence*>& db,
		size_t k,
		size_t fragLength,
		const SubstitutionMatrix & matrix,
		uint maxResults,
		string& outFile //
	) {
		const uint Q = query.size();
		const uint R = db.size();
		ofstream out( outFile );

		auto getFragCount = [k, fragLength]( const Sequence* seq ) -> size_t {
			auto kmerCount = seq->Seq().size() + 1 - k;
			auto fragCount = Fragment::GetCount( kmerCount, fragLength );
			return fragCount;
		};

		size_t maxQueryFragCount = Util::Max( query.begin(), query.end(), getFragCount, 0 );
		size_t maxDbFragCount = Util::Max( db.begin(), db.end(), getFragCount, 0 );

#if USE_OMP
#pragma omp parallel
#endif
		{
			// cout << "Running with " << omp_get_num_threads() << " threads.\n";

			KnnVector<size_t, double> rankings( maxResults, -HUGE_VAL );
			vector<int> rowMinima( maxQueryFragCount );
			vector<int> colMinima( maxDbFragCount );

#if USE_OMP
#pragma omp for schedule(dynamic,1)
#endif
			for ( uint q = 0; q < Q; q++ ) {
				rankings.clear();

				const auto& querySeq = query[q]->Seq();

				if (querySeq.size() < k) {
					out << query[q]->IdString();
					out << " ___eol___ -100000\n";
					continue;
				}

				size_t m = querySeq.size() + 1 - k;
				size_t queryFragCount = Fragment::GetCount( m, fragLength );
				double queryStepSize = Fragment::GetRealStepSize( m, fragLength, queryFragCount );

				for ( uint r = 0; r < R; r++ ) {
					const auto& refSeq = db[r]->Seq();

					if (refSeq.size() < k) continue;

					size_t n = refSeq.size() + 1 - k;
#if false
#pragma omp critical
					{
						cout << "Processing reference sequence " << db[r]->IdString() << "\n";
						cout.flush();
					}
#endif
					size_t refFragCount = Fragment::GetCount( n, fragLength );
					double refStepSize = Fragment::GetRealStepSize( n, fragLength, refFragCount );


					std::fill_n( rowMinima.begin(), queryFragCount, numeric_limits<int>::max() );
					std::fill_n( colMinima.begin(), refFragCount, numeric_limits<int>::max() );

					auto update = [&rowMinima, &colMinima, queryStepSize, refStepSize]( size_t i, size_t j, int d ) {
						size_t i_ = (size_t) ( i / queryStepSize );
						size_t j_ = (size_t) ( j / refStepSize );

						if ( rowMinima[i_] > d ) rowMinima[i_] = d;
						if ( colMinima[j_] > d ) colMinima[j_] = d;
					};

					Simproj::ComputeKmerDistances(querySeq, refSeq, k, matrix, update);

					double distance = Agg( rowMinima, m, colMinima, n );

					// cout << query.at(q)->IdStr() << "\t" << seq.IdStr() << "\t" << distance << "\n";

					if ( rankings.canPush( distance ) ) {
						rankings.push( r, distance );
					}
				}

				rankings.sort();

#if USE_OMP
#pragma omp critical
#endif
				{
					out << query[q]->IdString();

					for ( auto& ranking : rankings ) {
						out << " " << db[ranking.second]->IdString() << " " << ( -ranking.first );
					}

					out << " ___eol___ -100000\n";
				}
			}

		}
	}
};

int main( int argc, char* argv[] ) {
	try {
		Args args( argc, argv );

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AASP::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << ( end_time - start_time ) << "s" << endl;

		return retCode;
	}
	catch ( Exception& ex ) {
		cout << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}
