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

struct AAClustGetRandomPrototypes {

	using DistanceFunction = KmerDistanceCache2;
	using Fasta = FastaSequence;
	using Seq = EncodedFastaSequence;
	using Proto = KmerClusterPrototype;
	using SelectionMode = enum { random, minhash };

public:
	static int Run() {
		Params p;

		cerr << arguments->ProgName() << " \\\n" << p;

		BlosumDifferenceFunction rawDistanceFunction(p.matrix);
		DistanceFunction distanceFunction(p.alphabet, &rawDistanceFunction);

		omp_set_num_threads(p.numThreads);

		auto dbSeqs = Load::Fasta(p.fastaFile, p.idIndex, p.alphabet);
		auto db = Load::Encoded( dbSeqs, -1, p.alphabet, p.wordLength,
			distanceFunction.CharsPerWord(), p.alphabet->DefaultSymbol());

		cerr << arguments->ProgName() << ": " << db.size() << " sequences loaded.\n";

		KmerIndex kmerIndex(db, p.wordLength);
		auto & kmers = kmerIndex.GetKmers();

		cerr << arguments->ProgName() << ": " << kmers.size() << " distinct k-mers indexed.\n";

		vector<FastaSequence *> protoSeqs;

		if (p.selectMode == SelectionMode::random) {
			size_t N = kmers.size();

			UniformIntRandom<size_t> rand(p.seed, 0, N - 1);

			ofstream protoStream(p.protoFile);

			for (uint i = 0; i < kmers.size() && i < p.numProtos; i++) {
				size_t next = i + rand() % (N - i);
				std::swap(kmers[i], kmers[next]);
				auto protoSeq = new FastaSequence("", kmers[i]->Word(), 0, p.alphabet);
				protoSeqs.push_back( protoSeq );
				Proto prototype(protoSeq, p.alphabet, p.wordLength, distanceFunction.CharsPerWord(), p.alphabet->DefaultSymbol());
				protoStream << &prototype;
			}
		}

		else if (p.selectMode == SelectionMode::minhash) {
			KnnVector<Kmer *, size_t> knn(p.numProtos, numeric_limits<Distance>::max());
			Substring::Hash hash;

#if USE_OMP
#pragma omp parallel
#endif
			{
				KnnVector<Kmer *, size_t> knnPerThread(p.numProtos, numeric_limits<Distance>::max());

#if USE_OMP
#pragma omp for
#endif
				for (uint i = 0; i < kmers.size(); i++) {
					auto hashCode = hash(kmers[i]->Substr());

					if (knnPerThread.canPush(hashCode)) {
						knn.push(kmers[i], hashCode);
					}
				}

#if USE_OMP
#pragma omp critical
#endif
				{
					for (auto pair : knnPerThread) {
						if (knn.canPush(pair.first)) {
							knn.push(pair.second, pair.first);
						}
					}
				}
			}

			ofstream protoStream(p.protoFile);

			for (auto x : knn) {
				auto protoSeq = new FastaSequence("", x.second->Word(), 0, p.alphabet);
				protoSeqs.push_back(protoSeq);
				Proto prototype(protoSeq, p.alphabet, p.wordLength, distanceFunction.CharsPerWord(), p.alphabet->DefaultSymbol());
				protoStream << &prototype;
			}

			Util::Free(db);
			Util::Free(dbSeqs);
			Util::Free(protoSeqs);
		}

		return 0;
	}

	struct Params {
		bool ok = true;
		string fastaFile;
		string protoFile;
		size_t numThreads = 7;
		size_t wordLength;
		size_t idIndex = 0;
		int seed = time(nullptr);
		size_t numProtos;
		SelectionMode selectMode;
		SimilarityMatrix * matrix;
		Alphabet * alphabet;

		Params() {
#define REQ(x,h) arguments->Required(x, #x, h)
#define OPT(x,h) arguments->Optional(x, #x, h)
#define STR(x) STR2(x)
#define STR2(x) #x

			REQ(fastaFile,
				"The name of a file that will be parsed as a FASTA file which contains a sequence database.");

			REQ(protoFile,
				"The name of a file that will be populated with a FASTA-formatted list of prototypes.");

			REQ(idIndex, "The 0-origin position of the sequence ID field in the pipe-separated definition line.");

			REQ(seed, "The seed for the random number generator.");

			REQ(numProtos, "The number of prototypes to select from the .");

			OPT(numThreads, "default value = " STR(DEFAULT_THREADS) ". The number of OpenMP threads to use in parallel regions.");

			REQ(wordLength, "The word length used for kmer tiling.");

			arguments->Required(alphabet, matrix);

			vector<string> selectVals{ "random", "minhash" };
			string selectMode;
			bool selectOk = false;

			REQ(selectMode, "The selection mode. (random|minhash).");
			String::ToLowerInPlace(selectMode);

			for (uint i = 0; i < selectVals.size(); i++) {
				if (selectMode == selectVals[i]) {
					selectOk = true;
					this->selectMode = SelectionMode(i);
					break;
				}
			}

			if (!selectOk) {
				cerr << "Bad selection mode.\n";
				arguments->Fail();
			}

			if (!arguments->Ok()) {
				arguments->Help();
				throw Exception("Error processing command line arguments.", FileAndLine);
			}
		}

		friend ostream & operator << (ostream & str, const Params & p) {
			str << "--fastaFile '" << p.fastaFile << "' \\\n"
				<< "--protoFile '" << p.protoFile << "' \\\n"
				<< "--numThreads '" << p.numThreads << "' \\\n"
				<< "--wordLength '" << p.wordLength << "' \\\n"
				<< "--idIndex '" << p.idIndex << "' \\\n"
				<< "--seed '" << p.seed << "' \\\n"
				<< "--numProtos '" << p.numProtos << "' \\\n"
				<< "--alphabet '" << p.alphabet << "' \\\n"
				<< "--" << p.matrix->id << "\n";

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
		int retCode = AAClustGetRandomPrototypes::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception &ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}


