#include <AlphabetHelper.hpp>
#include <Args.hpp>
#include <Assert.hpp>
#include <Delegates.hpp>
#include <KmerCodebook.hpp>
#include <KmerDistanceCache.hpp>
#include <Edge.hpp>
#include <FastaSequence.hpp>
#include <Random.hpp>
#include <Kmer.hpp>
#include <KmerCluster.hpp>
#include <TestFramework.h>
#include <OmpTimer.h>
#include <FileUtil.hpp>
#include <kNearestNeighbours.hpp>

#include <bitset>
#include <cstdio>
#include <omp.h>
#include <set>
#include <DataLoader.hpp>

using namespace QutBio;
using namespace std;

// Address of singleton argument table.
Args *arguments;

struct GetRandomSubsetFasta {

	using DistanceFunction = KmerDistanceCache2;
	using Fasta = FastaSequence;
	using Seq = EncodedFastaSequence;
	using Proto = KmerClusterPrototype;
	using SelectionMode = enum { random, minhash };

public:
	static int Run() {
		Params p;

		cerr << arguments->ProgName() << " \\\n" << p;

		auto db = Load::Fasta(p.fastaFile, 0, Alphabets::DEFAULT());
		cerr << arguments->ProgName() << ": " << db.size() << " sequences loaded.\n";

		size_t N = db.size();

		UniformIntRandom<size_t> rand(p.seed, 0, N - 1);

		ofstream outStream(p.outFile);

		for (uint i = 0; i < db.size() && i < p.sampleSize; i++) {
			size_t next = i + rand() % (N - i);
			std::swap(db[i], db[next]);
			outStream << db[i];
		}

		return 0;
	}

	struct Params {
		bool ok = true;
		string fastaFile;
		string outFile;
		int seed = time(nullptr);
		size_t sampleSize;

		Params() {
#define REQ(x,h) arguments->Required(x, #x, h)
#define OPT(x,h) arguments->Optional(x, #x, h)
#define STR(x) STR2(x)
#define STR2(x) #x

			REQ(fastaFile,
				"The name of a file that will be parsed as a FASTA file which contains a sequence database.");

			REQ(outFile,
				"The name of a file that will be populated with a FASTA-formatted list of prototypes.");

			REQ(seed, "The seed for the random number generator.");

			REQ(sampleSize, "The number of prototypes to select from the .");

			if (!arguments->Ok()) {
				arguments->Help();
				throw Exception("Error processing command line arguments.", FileAndLine);
			}
		}

		friend ostream & operator << (ostream & str, const Params & p) {
			str << "--fastaFile '" << p.fastaFile << "' \\\n"
				<< "--outFile '" << p.outFile << "' \\\n"
				<< "--seed '" << p.seed << "' \\\n"
				<< "--sampleSize '" << p.sampleSize << "\n";

			return str;
		}
	};
};

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

int main(int argc, char *argv[]) {
	try {
		Args args(argc, argv);

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = GetRandomSubsetFasta::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception &ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}


