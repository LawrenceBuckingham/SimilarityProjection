#include "AlphabetHelper.hpp"
#include "Args.hpp"
#include "Assert.hpp"
#include "Delegates.hpp"
#include "KmerCodebook.hpp"
#include "KmerDistanceCache.hpp"
#include "Edge.hpp"
#include "FastaSequence.hpp"
#include "Random.hpp"
#include "Kmer.hpp"
#include "OmpTimer.h"
#include "SparseSignature.hpp"
#include "DataLoader.hpp"

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
	using Seq = EncodedFastaSequence;
	using Fasta = FastaSequence;
	using Proto = KmerClusterPrototype;

	struct Counter {
		Proto *proto;
		size_t freq;

		void increment() {
#pragma omp atomic
			freq++;
		}

		// Reverse-order comparator because I want the results to appear in descending order.
		bool operator<(const Counter & other) {
			return this->freq > other.freq;
		}
	};

	static int Run() {
		Params parms;

		if (!parms.ok) {
			return 1;
		}

		if (parms.numThreads > 0) {
			omp_set_num_threads(parms.numThreads);
		}

		auto protoSeqs = Load::Fasta(parms.protoFile, 0, Alphabets::AA());
		auto protos = Load::Prototypes(protoSeqs, Alphabets::AA(), parms.wordLength, 2);

		auto dbSigs = Signature::Read(parms.dbSigs);

		OMP_TIMER_DECLARE(rank);
		OMP_TIMER_START(rank);

		vector<Counter> histogram(protos.size());

#pragma omp parallel for
		for (uint i = 0; i < protos.size(); i++) {
			histogram[i].proto = protos[i];
		}

		const int N = dbSigs.size();

#pragma omp parallel for
		for (int i = 0; i < N; i++) {
			auto sig = dbSigs[i];

			for (auto index : *sig) {
				histogram[index].increment();
			}
		}

		std::sort(histogram.begin(), histogram.end());

		ofstream out(parms.outFile);

		for (auto & i : histogram) {
			out << i.freq << "\t" << i.proto->IdStr() << "\n";
		}

		OMP_TIMER_END(rank);

		Util::Free(protos);
		Util::Free(protoSeqs);

		return 0;
	}


	struct Params {
	public:
		string dbSigs;
		string protoFile;
		string outFile;
		size_t numThreads = -1;
		size_t wordLength;
		bool ok = true;

		Params() {

			if (arguments->IsDefined("help")) {
				vector<string> text{
"AAClustCountIndices: Emits a histogram showing the number of references to",
"             each prototype in a signature database. The entries are sorted",
"             in descending frequency order.",
"",
"",
"--help       Gets this text.",
"",
"--dbSigs     Required. The name of the file which contains signatures for the ",
"             reference sequences.",
"",
"--protoFile  Required. The name of the FASTA-formatted file which contains ",
"             prototype kmer definitions.",
"",
"--outFile    Required. The name of the output file. This will be a tab ",
"             delimited document with records containing two fields: the prototype ",
"             sequence ID and frequency with which that prototype appears in a ",
"             signature.",
"",
"--wordLength Required. The k-mer length.",
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
					"; running with default value.\n";
			}

			if (!arguments->Get("dbSigs", dbSigs)) {
				cerr << arguments->ProgName() << ": error - required argument '--dbSigs' not set.\n";
				ok = false;
			}

			if (!arguments->Get("wordLength", wordLength)) {
				cerr << arguments->ProgName() << ": error - required argument '--wordLength' not set.\n";
				ok = false;
			}

			if (!arguments->Get("protoFile", protoFile)) {
				cerr << arguments->ProgName() << ": error - required argument '--protoFile' not set.\n";
				ok = false;
			}

			if (!arguments->Get("outFile", outFile)) {
				cerr << arguments->ProgName() << ": note - required argument '--outFile' not provided.\n";
				ok = false;
			}

			if (outFile == dbSigs || outFile == protoFile) {
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
