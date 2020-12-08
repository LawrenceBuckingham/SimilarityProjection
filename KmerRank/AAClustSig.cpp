#include "AlphabetHelper.hpp"
#include "Args.hpp"
#include "Assert.hpp"
#include "DataLoader.hpp"
#include "Delegates.hpp"
#include "KmerCodebook.hpp"
#include "KmerDistanceCache.hpp"
#include "Edge.hpp"
#include "FastaSequence.hpp"
#include "Random.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"
#include "TestFramework.h"
#include "OmpTimer.h"
#include "BitSet.hpp"
#include "kNearestNeighbours.hpp"
#include "Ranking.hpp"
#include "SparseSignature.hpp"
#include <cstdio>
#include <stdlib.h>
#include <omp.h>
#include <set>
#include <vector>
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <utmpx.h>

using namespace QutBio;
using namespace std;

// Singletons.
Args *arguments;
mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;
function<const FastaSequence *(const string & seqId)> SparseSignature::Lookup;

struct AAClustSig {
public:
	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using pCluster = Cluster * ;
	using Signature = SparseSignature;

	static int Run() {
		Params parms;

		if (!parms.ok) {
			return 1;
		}

		if (parms.numThreads > 0) {
			omp_set_num_threads(parms.numThreads);
		}

		auto querySigs = Signature::Read(parms.querySigs);
		auto dbSigs = Signature::Read(parms.dbSigs);
		vector<vector<uint>> dbIndex;
		Signature::CreatePostingList(dbSigs, dbIndex);

		cerr << querySigs.size() << " query signatures loaded from '" << parms.querySigs << "'\n";
		cerr << dbSigs.size() << " query signatures loaded from '" << parms.dbSigs << "'\n";

		//for (auto & qs : querySigs) {
		//	cerr << "qs.indices.size() = " << qs->indices.size() << "\n";
		//}

		OMP_TIMER_DECLARE(rank);
		OMP_TIMER_START(rank);
		RankMerge(querySigs, dbSigs, dbIndex, parms.maxResults, parms.outFile);
		OMP_TIMER_END(rank);

		for (auto sig : querySigs) delete sig;
		for (auto sig : dbSigs) delete sig;

		return 0;
	}

	static void RankMerge(
		const vector<Signature *> &queries,
		const vector<Signature *> &database,
		const vector<vector<uint>> &dbIndex,
		uint maxResults,
		string &outFile //
	) {
		//cerr << "RankMerge\n";

		uint Q = queries.size();

#define INTERLEAVE 1

#if INTERLEAVE
		ofstream out(outFile);
#else
		KnnVector<size_t, double> exemplar(maxResults);
		vector<KnnVector<size_t, double>> allRankings(Q, exemplar);
#endif

#if USE_OMP
#pragma omp parallel
#endif
		{

#if INTERLEAVE
			KnnVector<size_t, double> rankings(maxResults, HUGE_VAL);
#endif
			BitSet processed(database.size());

#if USE_OMP
#pragma omp for
#endif
			for (uint q = 0; q < Q; q++) {
#if ! INTERLEAVE
				auto & rankings = allRankings[q];
#endif
				
				//cerr << "q = " << q << "\n";
				//cerr << "queryIndices.size() = " << queryIndices.size() << "\n";

				rankings.clear();
				processed.Clear();

				for (uint c : *queries[q]) {
					for (uint d : dbIndex[c]) {
						if (!processed.Contains(d)) {
							processed.Insert(d);
							double distance = 1.0 - queries[q]->Similarity(*database[d]);

							// cerr << queries[q]->id << " --> " << database[d]->id << " = " << distance << "\n";

							if (rankings.canPush(distance)) {
								rankings.push(d, distance);
							}
						}
					}
				}

				rankings.sort();

#if INTERLEAVE
#if USE_OMP
#pragma omp critical
#endif
				{
					out << queries[q]->Sequence()->IdStr();

					for (auto & ranking : rankings) {
						out << " " << database[ranking.second]->Sequence()->IdStr() << " " << (-ranking.first);
					}

					out << " ___eol___ -100000\n";
				}
#endif
			}

		}
#if ! INTERLEAVE
		{
			ofstream out(outFile);

			for (uint q = 0; q < Q; q++) {
				out << queries[q]->id;

				for (auto & ranking : allRankings[q]) {
					out << " " << database[ranking.second]->id << " " << (-ranking.first);
				}

				out << " ___eol___ -100000\n";
			}
		}
#endif
	}

	struct Params {
	public:
		string dbSigs;
		string querySigs;
		string outFile;
		size_t numThreads = 8;
		bool ok = true;
		uint maxResults = 1000;

		Params() {

			if (arguments->IsDefined("help")) {
				vector<string> text{
"AAClustSig:  Ranks the top K reference signatures for each sequence in a ",
"             query dataset. This program uses the Jaccard Similarity between ",
"             signatures of reference and query as a proxy for biological ",
"             sequence similarity. Requires precomputed signatures (use ",
"             AAClustSigEncode) but no fasta files, prototypes, or  clusters.",
"",
"--help       Gets this text.",
"",
"--dbSigs     Required. The name of the file which contains signatures for the ",
"             reference sequences.",
"",
"--querySigs  Required. The name of the file which contains signatures for the ",
"             query sequences.",
"",
"--outFile    Required. The name of the output file. This will be a CSV ",
"             document with records containing two fields: the prototype ",
"             sequence ID and information gain.",
"",
"--numThreads Optional; default value = '# cores'. The number of OpenMP ",
"             threads to use in parallel regions.",
"",
				};

				for (auto s : text) {
					cerr << s << "\n";
				}
			}

			if (!arguments->Get("numThreads", numThreads)) {
				cerr << arguments->ProgName() << ": note - optional argument '--numThreads' not set"
					"; running with default value "
					<< numThreads << ".\n";
			}

			if (!arguments->Get("dbSigs", dbSigs)) {
				cerr << arguments->ProgName() << ": error - required argument '--dbSigs' not set.\n";
				ok = false;
			}

			if (!arguments->Get("querySigs", querySigs)) {
				cerr << arguments->ProgName() << ": error - required argument '--querySigs' not set.\n";
				ok = false;
			}

			if (!arguments->Get("outFile", outFile)) {
				cerr << arguments->ProgName() << ": note - required argument '--outFile' not provided.\n";
				ok = false;
			}

			if (!arguments->Get("maxResults", maxResults)) {
				cerr << arguments->ProgName() << ": note - optional argument '--maxResults' not set"
					"; running with default value "
					<< maxResults << ".\n";
			}

			if (outFile == dbSigs || outFile == querySigs) {
				cerr << arguments->ProgName() << ": Output file " << outFile << " will overwrite one of your input files.\n";
				ok = false;
			}
		}
	};
};

int main(int argc, char *argv[]) {
	try {
		Args args(argc, argv);

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AAClustSig::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception &ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}
