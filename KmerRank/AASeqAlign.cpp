#include "AlphabetHelper.hpp"
#include "Args.hpp"
#include "Assert.hpp"
#include "Delegates.hpp"
#include "KmerCodebook.hpp"
#include "KmerDistanceCache.hpp"
#include "Edge.hpp"
#include "FastaSequence.hpp"
#include "HBRandom.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"
#include "TestFramework.h"
#include "OmpTimer.h"
#include "BitSet.hpp"
#include "kNearestNeighbours.hpp"
#include "Ranking.hpp"
#include <cstdio>
#include <stdlib.h>
#include <omp.h>
#include <set>
#include <vector>
#include "DataLoader.hpp"

#include "../seq-align/src/needleman_wunsch.h"
#include "../seq-align/src/smith_waterman.h"

using namespace QutBio;
using namespace std;

// Singletons.
Args *arguments;
mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

struct AASeqAlign {
public:

	struct Mode : public EnumBase {
	private:
		Mode(string literal, int value) : EnumBase(literal, value) {}

		static std::mutex & m() {
			static std::mutex mut;
			return mut;
		}

	public:
		static Mode * Needle() {
			std::unique_lock < mutex > lck{ m() };
			static Mode value("needle", 0);
			return &value;
		}

		static Mode * Water() {
			std::unique_lock < mutex > lck{ m() };
			static Mode value("water", 0);
			return &value;
		}

		static vector<EnumBase *> Values() {
			std::vector<EnumBase *> result{
				Needle(),
				Water()
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
		int idIndex = 0;
		bool ok = true;
		int maxResults = 500;
		Mode * mode = Mode::Needle();

		Params() {

			if (arguments->IsDefined("help")) {
				vector<string> text{
"AASeqAlign: Uses Needleman-Wunsch or Smith-Waterman to get the top maxResults ",
"               matches for a list of queries from a database.",
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
"--mode         Optional; default value = 'Needle'. The alignment algorithm to ",
"               apply. Valid values are: 'Needle', 'Water'. 'Needle' is ",
"               Needleman-Wunsch; 'Water' is Smith-Waterman. The present ",
"               version is hard-coded to use BLOSUM62, gap-open = 10; gap-extend ",
"               = 1.",
				};

				for (auto s : text) {
					if (s[0] == '-' && s[1] == '-') {
						cerr << "\n";
					}
					cerr << s << "\n";
				}
			}

			if (!arguments->Get("dbFile", this->dbFile)) {
				cerr << arguments->ProgName() << ": error - required argument '--dbFile' not supplied.\n";
				ok = false;
			}

			if (!arguments->Get("queryFile", this->queryFile)) {
				cerr << arguments->ProgName() << ": error - required argument '--queryFile' not "
					"supplied.\n";
				ok = false;
			}

			if (!arguments->Get("idIndex", idIndex)) {
				cerr << arguments->ProgName() << ": error - required argument '--idIndex' not "
					"supplied.\n";
				ok = false;
			}

			if (!arguments->Get("numThreads", numThreads)) {
				cerr << arguments->ProgName() << ": note - optional argument '--numThreads' not set"
					"; running with default value "
					<< numThreads << ".\n";
			}

			if (!arguments->Get("maxResults", maxResults)) {
				cerr << arguments->ProgName() << ": note - optional argument '--maxResults' not set"
					"; running with default value "
					<< maxResults << ".\n";
			}

			if (!arguments->Get("mode", Mode::Values(), mode)) {
				cerr << arguments->ProgName() << ": note - optional argument '--mode' not set"
					"; running with default value "
					<< (*mode) << ".\n";
			}

			if (!arguments->Get("outFile", outFile)) {
				cerr << arguments->ProgName() << ": note - required argument '--outFile' not provided.\n";
				ok = false;
			}
		}
	};

	static int Run() {
		Params parms;

		if (!parms.ok) {
			return 1;
		}

		omp_set_num_threads(parms.numThreads);

		auto db = Load::Fasta(parms.dbFile, parms.idIndex, Alphabets::AA());
		cerr << arguments->ProgName() << ": " << db.size() << " reference sequences loaded from " << parms.dbFile << ".\n";

		auto query = Load::Fasta(parms.queryFile, parms.idIndex, Alphabets::AA());
		cerr << arguments->ProgName() << ": " << query.size() << " query sequences loaded from " << parms.queryFile << ".\n";

		OMP_TIMER_DECLARE(rankTime);
		OMP_TIMER_START(rankTime);
		if (parms.mode == Mode::Needle()) {
			RankNeedle(query, db, parms.maxResults, parms.outFile);
		}
		else if (parms.mode == Mode::Water()) {
			RankWater(query, db, parms.maxResults, parms.outFile);
		}
		else {
			throw Exception("Whatever mode you entered is not implemented in the present version.", FileAndLine);
		}
		OMP_TIMER_END(rankTime);

		cerr << "Ranking completed in " << OMP_TIMER(rankTime) << "s.\n";

		return 0;
	}

	static void RankNeedle(
		vector<FastaSequence *> &query,
		vector<FastaSequence *> &db,
		uint maxResults,
		string &outFile //
	) {
		cerr << "Ranking sequences via Needleman-Wunsch\n";

		const uint Q = query.size();
		ofstream out(outFile);

#pragma omp parallel
		{
			KnnVector<size_t, double> rankings(maxResults, -HUGE_VAL);
			nw_aligner_t *nw = needleman_wunsch_new();
			alignment_t *result = alignment_create(256);
			scoring_t scoring;

			scoring_system_BLOSUM62(&scoring);
			scoring.gap_open = -10;
			scoring.gap_extend = -1;

#pragma omp for
			for (uint q = 0; q < Q; q++) {
				auto & querySeq = query[q];

				rankings.clear();

				for (uint r = 0; r < db.size(); r++) {
					auto & refSeq = db[r];
					double distance = -Needle(*querySeq, *refSeq, nw, result, &scoring);

					if (rankings.canPush(distance)) {
						rankings.push(r, distance);
					}
				}

				rankings.sort();

#pragma omp critical
				{
					out << query[q]->IdStr();

					for (auto & ranking : rankings) {
						out << " " << db[ranking.second]->IdStr() << " " << (-ranking.first);
					}

					out << " ___eol___ -100000\n";
				}
			}

			needleman_wunsch_free(nw);
			alignment_free(result);
		}
	}

	static void RankWater(
		vector<FastaSequence *> &query,
		vector<FastaSequence *> &db,
		uint maxResults,
		string &outFile //
	) {
		cerr << "Ranking sequences via Smith-Waterman\n";

		const uint Q = query.size();
		ofstream out(outFile);

#pragma omp parallel
		{
			KnnVector<size_t, double> rankings(maxResults, -HUGE_VAL);
			sw_aligner_t *sw = smith_waterman_new();
			alignment_t *result = alignment_create(256);
			scoring_t scoring;

			scoring_system_BLOSUM62(&scoring);
			scoring.gap_open = -10;
			scoring.gap_extend = -1;

#pragma omp for
			for (uint q = 0; q < Q; q++) {
				auto & querySeq = query[q];

				rankings.clear();

				for (uint r = 0; r < db.size(); r++) {
					auto & refSeq = db[r];
					double distance = -Water(*querySeq, *refSeq, sw, result, &scoring);

					if (rankings.canPush(distance)) {
						rankings.push(r, distance);
					}
				}

				rankings.sort();

#pragma omp critical
				{
					out << query[q]->IdStr();

					for (auto & ranking : rankings) {
						out << " " << db[ranking.second]->IdStr() << " " << (-ranking.first);
					}

					out << " ___eol___ -100000\n";
				}
			}

			smith_waterman_free(sw);
			alignment_free(result);
		}
	}

	static double Needle(
		FastaSequence & querySeq,
		FastaSequence & refSeq,
		nw_aligner_t *nw,
		alignment_t *alignmentResult,
		scoring_t* scoring
	) {
		needleman_wunsch_align(querySeq.CharData().c_str(), refSeq.CharData().c_str(), scoring, nw, alignmentResult);

		return alignmentResult->score;
	}

	static double Water(
		FastaSequence & querySeq,
		FastaSequence & refSeq,
		sw_aligner_t *sw,
		alignment_t *alignmentResult,
		scoring_t* scoring
	) {
		smith_waterman_align(querySeq.CharData().c_str(), refSeq.CharData().c_str(), scoring, sw);

		if (smith_waterman_fetch(sw, alignmentResult) > 0) {
			return alignmentResult->score;
		}
		else {
			return -10 - std::max(querySeq.Sequence().size(), refSeq.Sequence().size());
		}
	}

};

int main(int argc, char *argv[]) {
	try {
		Args args(argc, argv);

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AASeqAlign::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception & ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}
