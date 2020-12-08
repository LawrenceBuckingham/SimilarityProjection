#include "AlphabetHelper.hpp"
#include "Args.hpp"
#include "Assert.hpp"
#include "DataLoader.hpp"
#include "Delegates.hpp"
#include "KmerCodebook.hpp"
#include "KmerDistanceCache.hpp"
#include "Edge.hpp"
#include "EncodedFastaSequence.hpp"
#include "Random.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"
#include "TestFramework.h"
#include "OmpTimer.h"
#include "BitSet.hpp"
#include "kNearestNeighbours.hpp"
#include "Ranking.hpp"
#include "FragmentAggregationMode.hpp"
#include "Simproj.hpp"

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

struct AAD2 {
public:

	struct D2Mode : public EnumBase {
	private:
		D2Mode( string literal, int value ) : EnumBase( literal, value ) {}
		static std::mutex& m() {
			static std::mutex mut;
			return mut;
		}

	public:
		static D2Mode* D2() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "d2", 0 );
			return &value;
		}

		static D2Mode* E() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "e", 1 );
			return &value;
		}

		static D2Mode* E_norm() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "e_norm", 2 );
			return &value;
		}

#if WANT_D2_STAR
		static D2Mode* D2S() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "d2s", 3 );
			return &value;
		}

		static D2Mode* D2S_missing() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "d2s_missing", 4 );
			return &value;
		}

		static D2Mode* D2S_observed() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "d2s_observed", 5 );
			return &value;
		}

		static D2Mode* D2Star() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "d2star", 6 );
			return &value;
		}

		static D2Mode* D2Star_missing() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "d2star_missing", 7 );
			return &value;
		}

#endif

		static D2Mode* Cosine() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "cosine", 8 );
			return &value;
		}

		static D2Mode* Jaccard() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "jaccard", 9 );
			return &value;
		}

		static D2Mode* Min() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "min", 10 );
			return &value;
		}

		static D2Mode* MinNormMin() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "min_norm_min", 11 );
			return &value;
		}

		static D2Mode* MinNormMax() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "min_norm_max", 11 );
			return &value;
		}

		static D2Mode* MinNormAvg() {
			std::unique_lock < mutex > lck{ m() };
			static D2Mode value( "min_norm_avg", 12 );
			return &value;
		}

		static vector<EnumBase*> Values() {
			std::vector<EnumBase*> result{
				D2(),
#if WANT_D2_STAR
				D2Star(),
				D2Star_missing(),
				D2S(),
				D2S_missing(),
				D2S_observed(),
#endif
				Cosine(),
				Jaccard(),
				E(),
				E_norm(),
				Min(),
				MinNormMin(),
				MinNormAvg()
			};
			return result;
		}
	};

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
		D2Mode* d2Mode = D2Mode::D2();
		size_t fragLength = string::npos;
		FragmentAggregationMode* fragMode = FragmentAggregationMode::BestOfBest();

		Params() {

			if ( arguments->IsDefined( "help" ) ) {
				vector<string> text{
"AAD2: Uses D2 to get the top maxResults matches for a list of queries from a ",
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
"--d2Mode       Optional; default value = 'D2'. The D2 variant to apply. Valid ",
"               values are: D2, D2S, D2Star, D2S_missing. 'D2' is scalar ",
"               product of word count vectors; 'D2S' is the scalar product of ",
"               standardised word count vectors; and 'D2Star', an alternative ",
"               normalised form. Refer to the paper 'Alignment-Free Sequence ",
"               Comparison (I): Statistics and Power, by GESINE REINERT, DAVID ",
"               CHEW, FENGZHU SUN, and MICHAEL S.WATERMAN' for further ",
"               information. I have also implemented an 'optimised' form of ",
"               D2S copied from JD2Stat ",
"               (http://bioinformatics.org.au/tools/jD2Stat/), which overlooks ",
"               all k-mers which are not present in both sequences; this is ",
"               not the same statistic as Reinert et al implemented in their ",
"               reference code. However, if we allow that a k-mer might be ",
"               considered to be a missing observation, and then allow that ",
"               missing observations might be replaced by the expected value, ",
"               then I can see an interpretation which makes this correct. The ",
"               optimised version of D2S is called 'D2S_missing' in this ",
"               program. Similarly, the optimised version of D2* is called ",
"               'D2Star_missing'. In addition, I have tried an in-between ",
"               version called 'D2S_observed', which penalises words that are ",
"               present in one string but not the other by including them in ",
"               the calculation. As this nightmare unfolds, I continue to ",
"               explore ways of measuring the similarity between two bags of ",
"               words.",
"               Use 'cosine' to compute cosine similarity (normalise ",
"               histograms prior to D2).",
"               Use 'jaccard' to ignore term counts, and compute Jaccard index.",
"               Use 'e' to compute squared euclidean distance.",
"               Use 'e_norm' to compute squared euclidean distance of ",
"               normalised histograms.",
"--fragMode     Optional. Fragment aggregation mode: one of the following:",
"               BestOfBest - the overall minimum fragment distance.",
"               Hausdorff - the maximum of the minimum one-way fragment ",
"               distances.",
"               HausdorffAverage - the maximum of the average one-way minimum ",
"               fragment distances.",
"               HausedorffAverageAverage - the average of the average one-way ",
"               minimum fragment distances.",
"               Default value is BestOfBest.",
"--fragLength   Optional length of fragment. Default value is a very large ",
"               value, which effectively causes the sequence as a whole to be ",
"               taken as a single fragment. In practice, some fragments will be",
"               stretched by 1 to ensure that the entire sequence is covered ",
"               evenly.",
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
					"; running with default value "
					<< numThreads << ".\n";
			}

			if ( !arguments->Get( "maxResults", maxResults ) ) {
				cout << arguments->ProgName() << ": note - optional argument '--maxResults' not set"
					"; running with default value "
					<< maxResults << ".\n";
			}

			if ( !arguments->Get( "d2Mode", D2Mode::Values(), d2Mode ) ) {
				cout << arguments->ProgName() << ": note - optional argument '--d2Mode' not set"
					"; running with default value "
					<< ( *d2Mode ) << ".\n";
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
					<< "default value " << fragMode->Name() << " will be used.\n";
			}

			if ( !arguments->Get( "fragLength", fragLength ) ) {
				cout << arguments->ProgName() << ": note - optional argument '--fragLength' not provided.\n"
					<< "default value " << fragLength << " will be used.\n";
			}
		}
	};

	struct TermFreqRecord {
		using T = TermFreqRecord;

		size_t first;
		double second;

		TermFreqRecord( const Substring& key, double value ) : first( key.HashCode() ), second( value ) {}

		friend bool operator<( const T& lhs, const T& rhs ) {
			return lhs.first < rhs.first;
		}
	};

	struct TermFreqVector : public vector<TermFreqRecord> {
		double& operator[] ( const Substring& key ) {
			for ( auto& tuple : *this ) {
				if ( tuple.first == key.HashCode() ) {
					return tuple.second;
				}
			}

			emplace_back( key, 0.0 );
			return back().second;
		}

		friend
		ostream & operator<<( ostream & out, const TermFreqVector & v) {
			bool deja = false;

			for ( const auto& p : v ) {
				if (deja) out << ",";
				deja = true;
				out << "(" << p.first << "," << p.second << ")";
			}

			return out;
		}
	};

	static int Run() {
		Params p;

		if ( !p.ok ) {
			return 1;
		}

		auto alphabet = Alphabets::AA();

		omp_set_num_threads( p.numThreads );

		auto dbSeqs = Load::Fasta( p.dbFile, p.idIndex, alphabet );
		auto db = Load::Encoded( dbSeqs, -1, alphabet, p.wordLength, 1, alphabet->DefaultSymbol() );
		cout << arguments->ProgName() << ": " << db.size() << " reference sequences loaded.\n";

		auto querySeqs = Load::Fasta( p.queryFile, p.idIndex, alphabet );
		auto query = Load::Encoded( querySeqs, -1, alphabet, p.wordLength, 1, alphabet->DefaultSymbol() );
		cout << arguments->ProgName() << ": " << query.size() << " query sequences loaded.\n";

		auto symbolHistogram = FastaSequence::GetSymbolHistogram( dbSeqs );
		KmerHashIndex idx( db, p.wordLength, alphabet );
		cout << arguments->ProgName() << ": " << idx.size() << " k-mers indexed from database.\n";

		vector<vector<TermFreqVector>> dbTerms( db.size() );
		vector<TermFreqVector> dbSummaryTfv( db.size() );
		CreateTermVectors( db, dbTerms, dbSummaryTfv, p, alphabet );
		cout << arguments->ProgName() << ": Term Occurrence vectors generated from database.\n";

		vector<vector<TermFreqVector>> queryTerms( query.size() );
		vector<TermFreqVector> querySummaryTfv( query.size() );
		CreateTermVectors( query, queryTerms, querySummaryTfv, p, alphabet );
		cout << arguments->ProgName() << ": Term Occurrence vectors generated from query set.\n";

		//for ( auto & bag: queryBags ) {
		//	for ( auto &x: bag ) {
		//		cout << x.first << " --> " << x.second << "\n";
		//	}
		//	cout << "\n";
		//}

		OMP_TIMER_DECLARE( rankTime );
		OMP_TIMER_START( rankTime );

		if ( p.d2Mode == D2Mode::D2() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<D2, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<D2, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<D2, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<D2, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
#if WANT_D2_STAR
		else if ( p.d2Mode == D2Mode::D2S() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				RankP<D2S, BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				RankP<D2S, HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				RankP<D2S, HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				RankP<D2S, Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );
		}
#endif
		else if ( p.d2Mode == D2Mode::E() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<E, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<E, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<E, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<E, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::E_norm() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<E, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<E, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<E, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<E, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
#if WANT_D2_STAR
		else if ( p.d2Mode == D2Mode::D2S_missing() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				RankP<D2S_missing, BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				RankP<D2S_missing, HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				RankP<D2S_missing, HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				RankP<D2S_missing, Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::D2Star() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				RankP<D2Star, BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				RankP<D2Star, HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				RankP<D2Star, HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				RankP<D2Star, Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::D2Star_missing() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				RankP<D2Star_missing, BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				RankP<D2Star_missing, HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				RankP<D2Star_missing, HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				RankP<D2Star_missing, Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::D2S_observed() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				RankP<D2S_observed, BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				RankP<D2S_observed, HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				RankP<D2S_observed, HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				RankP<D2S_observed, Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, symbolHistogram, p.outFile );
		}
#endif
		else if ( p.d2Mode == D2Mode::Cosine() ) {
			// normalisation has already been taken care of.
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<D2, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<D2, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<D2, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<D2, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::Jaccard() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<Jaccard, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<Jaccard, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<Jaccard, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<Jaccard, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::Min() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<Min, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<Min, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<Min, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<Min, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::MinNormMin() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<MinNormMin, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<MinNormMin, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<MinNormMin, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<MinNormMin, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::MinNormMax() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<MinNormMax, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<MinNormMax, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<MinNormMax, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<MinNormMax, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
		else if ( p.d2Mode == D2Mode::MinNormAvg() ) {
			if ( p.fragMode == FragmentAggregationMode::BestOfBest() )
				Rank<MinNormAvg, Simproj::BestOfBest>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverageAverage() )
				Rank<MinNormAvg, Simproj::HausdorffAverageAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::HausdorffAverage() )
				Rank<MinNormAvg, Simproj::HausdorffAverage>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );

			else if ( p.fragMode == FragmentAggregationMode::Hausdorff() )
				Rank<MinNormAvg, Simproj::Hausdorff>( query, queryTerms, querySummaryTfv, db, dbTerms, dbSummaryTfv, idx, p.maxResults, p.outFile );
		}
		else {
			throw Exception( "Whatever d2 mode you entered is not implemented in the present version.", FileAndLine );
		}
		OMP_TIMER_END( rankTime );

#if USE_OMP
		cout << "Ranking completed in " << OMP_TIMER( rankTime ) << "s.\n";
#endif

		Util::Free( query );
		Util::Free( querySeqs );
		Util::Free( db );
		Util::Free( dbSeqs );

		return 0;
	}

	static void CreateTermVectors(
		vector<EncodedFastaSequence*>& db,
		vector<vector<TermFreqVector>>& terms,
		vector<TermFreqVector>& summaryTerms,
		AAD2::Params& parms,
		const Alphabet * alphabet
	) {
		uint i = 0;

		for ( auto seq_ : db ) {
			auto& seq = *seq_;
			seq.position = i;

			TermFreqVector& summaryTfv = summaryTerms[i];
			vector<TermFreqVector>& frags = terms[i];

			i++;

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

			seq.SelectKmers( parms.wordLength, process );

			if ( parms.d2Mode == D2Mode::Cosine() || parms.d2Mode == D2Mode::E_norm() ) {
				for ( auto& tfv : frags ) {
					Normalise( tfv );
				}

				Normalise( summaryTfv );
			}

			for ( auto& tfv : frags ) {
				std::sort( tfv.begin(), tfv.end() );
			}

			std::sort( summaryTfv.begin(), summaryTfv.end() );
		}
	}

	static void Normalise( TermFreqVector& bag ) {
		double sumSquared = 0;

		for ( auto& term : bag ) {
			sumSquared += term.second * term.second;
		}

		double norm = sqrt( sumSquared );


		for ( auto& term : bag ) {
			term.second /= norm;
		}
	}

	using BagSimilarity = double( * )( const TermFreqVector& x, const TermFreqVector& y );

	using Aggregator = double( * )(
		const vector<double>& rowMinima, size_t rowCount,
		const vector<double>& colMinima, size_t colCount
		);

	using BagSimilarityProb = double( * )( const TermFreqVector& x, const TermFreqVector& y, const double p[256] );

	template<BagSimilarity cmp, Aggregator agg>
	static void Rank(
		vector<EncodedFastaSequence*>& query,
		vector<vector<TermFreqVector>>& queryBags,
		vector<TermFreqVector>& querySummaryTfvs,

		vector<EncodedFastaSequence*>& db,
		vector<vector<TermFreqVector>>& dbBags,
		vector<TermFreqVector>& dbSummaryTfvs,

		KmerHashIndex& dbIndex,
		uint maxResults,
		string& outFile //
	) {
		const uint Q = query.size();
		ofstream out( outFile );

		size_t maxQueryFragCount = GetMaxFragCount( queryBags );
		size_t maxDbFragCount = GetMaxFragCount( dbBags );

#if USE_OMP
#pragma omp parallel
#endif
		{
			// cout << "Running with " << omp_get_num_threads() << " threads.\n";

			KnnVector<size_t, double> rankings( maxResults, -HUGE_VAL );
			BitSet processed( db.size() );
			vector<double> rowMinima( maxQueryFragCount );
			vector<double> colMinima( maxDbFragCount );

#if USE_OMP
#pragma omp for schedule(dynamic,1)
#endif
			for ( uint q = 0; q < Q; q++ ) {
				auto& querySummaryTfv = querySummaryTfvs[q];

				// cout << querySummaryTfv << "\n";

				rankings.clear();
				processed.Clear();

				for ( auto c : querySummaryTfv ) {
					auto dbKmer = dbIndex[c.first];

					if ( dbKmer == 0 ) continue;

					// cout << dbKmer->Substr() << "\n";

					for ( auto instance : dbKmer->Instances() ) {
						auto& seq = instance.sequence;
						size_t d = seq.position;

						if ( !processed.Contains( d ) ) {
							processed.Insert( d );

							vector<TermFreqVector>& queryFrags = queryBags[q];
							vector<TermFreqVector>& dbFrags = dbBags[d];

							size_t m = queryFrags.size();
							size_t n = dbFrags.size();

							//cout << "m = " << m << endl;
							//cout << "n = " << n << endl;

							std::fill_n( rowMinima.begin(), m, numeric_limits<double>::max() );
							std::fill_n( colMinima.begin(), n, numeric_limits<double>::max() );

							for ( size_t i = 0; i < m; i++ ) {
								for ( size_t j = 0; j < n; j++ ) {
									double fragDist = cmp( queryFrags[i], dbFrags[j] );
									if ( rowMinima[i] > fragDist ) rowMinima[i] = fragDist;
									if ( colMinima[j] > fragDist ) colMinima[j] = fragDist;
								}
							}

							double distance = agg( rowMinima, m, colMinima, n );

							// cout << query.at(q)->IdStr() << "\t" << seq.IdStr() << "\t" << distance << "\n";

							if ( rankings.canPush( distance ) ) {
								rankings.push( d, distance );
							}
						}
					}
				}

				rankings.sort();

#if USE_OMP
#pragma omp critical
#endif
				{
					out << query[q]->IdStr();

					for ( auto& ranking : rankings ) {
						out << " " << db[ranking.second]->IdStr() << " " << ( -ranking.first );
					}

					out << " ___eol___ -100000\n";
				}
			}

		}
	}

	static size_t GetMaxFragCount( std::vector<std::vector<AAD2::TermFreqVector>>& queryBags ) {
		size_t maxQueryFragCount = 0;

		for ( auto& queryFragList : queryBags ) {
			size_t queryFragCount = queryFragList.size();

			if ( queryFragCount > maxQueryFragCount ) {
				maxQueryFragCount = queryFragCount;
			}
		}

		return maxQueryFragCount;
	}

	template<BagSimilarityProb cmp, Aggregator agg>
	static void RankP(
		vector<EncodedFastaSequence*>& query,
		vector<vector<TermFreqVector>>& queryBags,
		vector<TermFreqVector>& querySummaryTfvs,

		vector<EncodedFastaSequence*>& db,
		vector<vector<TermFreqVector>>& dbBags,
		vector<TermFreqVector>& dbSummaryTfvs,

		KmerHashIndex& dbIndex,
		uint maxResults,
		Histogram<byte>& symbolHist,
		string& outFile //
	) {
		// Convert histogram to efficient lookup table.
		double p[256] = { 0 };

		symbolHist.Normalise();
		auto symbols = symbolHist.GetKeys();

		for ( auto symbol : symbols ) {
			p[(size_t) symbol] = symbolHist[symbol];
		}

		const uint Q = query.size();
		ofstream out( outFile );

		size_t maxQueryFragCount = GetMaxFragCount( queryBags );
		size_t maxDbFragCount = GetMaxFragCount( dbBags );

#pragma omp parallel
		{
			KnnVector<size_t, double> rankings( maxResults, -HUGE_VAL );
			BitSet processed( db.size() );
			vector<double> rowMinima( maxQueryFragCount );
			vector<double> colMinima( maxDbFragCount );

#pragma omp for
			for ( uint q = 0; q < Q; q++ ) {
				auto& querySummaryTfv = querySummaryTfvs[q];

				rankings.clear();
				processed.Clear();

				for ( auto c : querySummaryTfv ) {
					auto dbKmer = dbIndex[c.first];

					if ( dbKmer == 0 ) continue;

					for ( auto instance : dbKmer->Instances() ) {
						auto& seq = instance.sequence;
						size_t d = seq.position;

						if ( !processed.Contains( d ) ) {
							processed.Insert( d );

							vector<TermFreqVector>& queryFrags = queryBags[q];
							vector<TermFreqVector>& dbFrags = dbBags[q];

							size_t m = queryFrags.size();
							size_t n = dbFrags.size();

							std::fill_n( rowMinima.begin(), m, numeric_limits<double>::max() );
							std::fill_n( colMinima.begin(), n, numeric_limits<double>::max() );

							for ( size_t i = 0; i < m; i++ ) {
								for ( size_t j = 0; j < n; j++ ) {
									double fragDist = cmp( queryFrags[i], dbFrags[j], p );
									if ( rowMinima[i] > fragDist ) rowMinima[i] = fragDist;
									if ( colMinima[j] > fragDist ) colMinima[j] = fragDist;
								}
							}

							double distance = agg( rowMinima, m, colMinima, n );
							if ( rankings.canPush( distance ) ) {
								rankings.push( d, distance );
							}
						}
					}
				}

				rankings.sort();

#pragma omp critical
				{
					out << query[q]->IdStr();

					for ( auto& ranking : rankings ) {
						out << " " << db[ranking.second]->IdStr() << " " << ( -ranking.first );
					}

					out << " ___eol___ -100000\n";
				}
			}
		}
	}

	static double D2( const TermFreqVector& a, const TermFreqVector& b ) {
		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double sum = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				sum += i->second * j->second;
				i++;
				j++;
			}
		}

		return -sum;
	}

	static double E( const TermFreqVector& a, const TermFreqVector& b ) {
		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double sum = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				sum += i->second * i->second;
				i++;
			}
			else if ( y < x ) {
				sum += j->second * j->second;
				j++;
			}
			else {
				double t = ( i->second - j->second );
				sum += t * t;
				i++;
				j++;
			}
		}

		return sum;
	}

	static double Min( const TermFreqVector& a, const TermFreqVector& b ) {
		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double score = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				score += std::min( i->second, j->second );
				i++;
				j++;
			}
		}

		return -score;
	}

	static double MinNormMin( const TermFreqVector& a, const TermFreqVector& b ) {
		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double score = 0, aLen = 0, bLen = 0;

		for ( auto p = i; p != m; p++ ) {
			aLen += p->second;
		}

		for ( auto p = j; p != j; p++ ) {
			bLen += p->second;
		}

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				score += std::min( i->second, j->second );
				i++;
				j++;
			}
		}

		auto len = std::min( aLen, bLen );

		return len == 0 ? 0 : -score / len;
	}

	static double MinNormMax( const TermFreqVector& a, const TermFreqVector& b ) {
		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double score = 0, aLen = 0, bLen = 0;

		for ( auto p = i; p != m; p++ ) {
			aLen += p->second;
		}

		for ( auto p = j; p != j; p++ ) {
			bLen += p->second;
		}

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				score += std::min( i->second, j->second );
				i++;
				j++;
			}
		}

		auto len = std::max( aLen, bLen );

		return len == 0 ? 0 : -score / len;
	}

	static double MinNormAvg( const TermFreqVector& a, const TermFreqVector& b ) {
		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double score = 0, aLen = 0, bLen = 0;

		for ( auto p = i; p != m; p++ ) {
			aLen += p->second;
		}

		for ( auto p = j; p != j; p++ ) {
			bLen += p->second;
		}

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				score += std::min( i->second, j->second );
				i++;
				j++;
			}
		}

		auto len = ( aLen + bLen ) / 2;

		return len == 0 ? 0 : -score / len;
	}

	static double Jaccard( const TermFreqVector& a, const TermFreqVector& b ) {
		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double union_ = 0, intersect = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			union_++;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				i++;
				j++;
				intersect++;
			}
		}

		return 1 - intersect / union_;
	}

#if WANT_D2_STAR

	static double D2S( const TermFreqVector& a, const TermFreqVector& b, const double p[256] ) {
		double M = Sum( a );
		double N = Sum( b );

		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double sumScores = 0;
		double sumProbWordsSeen = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				double prob = Probability( i->first, p );
				sumProbWordsSeen += prob;
				double Xw = i->second - M * prob;
				double Yw = -N * prob;
				sumScores += Xw * Yw / sqrt( Xw * Xw + Yw * Yw );
				i++;
			}
			else if ( y < x ) {
				double prob = Probability( j->first, p );
				sumProbWordsSeen += prob;
				double Xw = -M * prob;
				double Yw = j->second - N * prob;
				sumScores += Xw * Yw / sqrt( Xw * Xw + Yw * Yw );
				j++;
			}
			else {
				double prob = Probability( i->first, p );
				sumProbWordsSeen += prob;
				double Xw = i->second - M * prob;
				double Yw = j->second - N * prob;
				sumScores += Xw * Yw / sqrt( Xw * Xw + Yw * Yw );
				i++;
				j++;
			}
		}

		return -( sumScores + M * N * ( 1 - sumProbWordsSeen ) / sqrt( M * M + N * N ) );
	}

	static double D2S_observed( const TermFreqVector& a, const TermFreqVector& b, const double p[256] ) {
		double M = Sum( a );
		double N = Sum( b );

		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double sumScores = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				double prob = Probability( i->first, p );
				double Xw = i->second - M * prob;
				double Yw = -N * prob;
				sumScores += Xw * Yw / sqrt( Xw * Xw + Yw * Yw );
				i++;
			}
			else if ( y < x ) {
				double prob = Probability( j->first, p );
				double Xw = -M * prob;
				double Yw = j->second - N * prob;
				sumScores += Xw * Yw / sqrt( Xw * Xw + Yw * Yw );
				j++;
			}
			else {
				double prob = Probability( i->first, p );
				double Xw = i->second - M * prob;
				double Yw = j->second - N * prob;
				sumScores += Xw * Yw / sqrt( Xw * Xw + Yw * Yw );
				i++;
				j++;
			}
		}

		return -sumScores;
	}

	static double D2S_missing( const TermFreqVector& a, const TermFreqVector& b, const double p[256] ) {
		double M = Sum( a );
		double N = Sum( b );

		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double sumScores = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				double prob = Probability( i->first, p );
				double Xw = i->second - M * prob;
				double Yw = j->second - N * prob;
				sumScores += Xw * Yw / sqrt( Xw * Xw + Yw * Yw );
				i++;
				j++;
			}
		}

		return -sumScores;
	}

	static double D2Star( const TermFreqVector& a, const TermFreqVector& b, const double p[256] ) {
		double M = Sum( a );
		double N = Sum( b );

		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double sumScores = 0;
		double sumProbWordsSeen = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				double prob = Probability( i->first, p );
				sumProbWordsSeen += prob;
				double Xw = i->second - M * prob;
				double Yw = -N * prob;
				sumScores += Xw * Yw / ( prob * sqrt( M * N ) );
				i++;
			}
			else if ( y < x ) {
				double prob = Probability( j->first, p );
				sumProbWordsSeen += prob;
				double Xw = -M * prob;
				double Yw = j->second - N * prob;
				sumScores += Xw * Yw / ( prob * sqrt( M * N ) );
				j++;
			}
			else {
				double prob = Probability( i->first, p );
				sumProbWordsSeen += prob;
				double Xw = i->second - M * prob;
				double Yw = j->second - N * prob;
				sumScores += Xw * Yw / ( prob * sqrt( M * N ) );
				i++;
				j++;
			}
		}

		return -( sumScores + sqrt( M * N ) * ( 1 - sumProbWordsSeen ) );
	}

	static double D2Star_missing( const TermFreqVector& a, const TermFreqVector& b, const double p[256] ) {
		double M = Sum( a );
		double N = Sum( b );

		auto i = a.begin();
		auto m = a.end();

		auto j = b.begin();
		auto n = b.end();

		double sumScores = 0;

		while ( i != m && j != n ) {
			auto x = i->first;
			auto y = j->first;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				double prob = Probability( i->first, p );
				double Xw = i->second - M * prob;
				double Yw = j->second - N * prob;
				sumScores += Xw * Yw / ( prob * sqrt( M * N ) );
				i++;
				j++;
			}
		}

		return -sumScores;
	}

	static double Probability( const Substring& substr, const double p[256] ) {
		double prob = 1.0;
		auto chars = substr.Chars();

		for ( size_t i = 0; i < substr.Length(); i++ ) {
			prob *= p[chars[i]];
		}

		return prob;
	}

#endif

	static double Sum( const TermFreqVector& a ) {
		double sum = 0;

		for ( auto& i : a ) {
			sum += i.second;
		}

		return sum;
	}
};

int main( int argc, char* argv[] ) {
	try {
		Args args( argc, argv );

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AAD2::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << ( end_time - start_time ) << "s" << endl;

		return retCode;
	}
	catch ( Exception& ex ) {
		cout << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}
