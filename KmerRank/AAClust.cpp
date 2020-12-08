#include "Args.hpp"
#include "Assert.hpp"
#include "Delegates.hpp"
#include "SubstitutionMatrix.hpp"
#include "Edge.hpp"
#include "Random.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"
#include "TestFramework.h"
#include "OmpTimer.h"
#include "KmerClusterPrototype.hpp"
#include "FileUtil.hpp"
#include "DataLoader.hpp"

#include <bitset>
#include <cstdio>
#include <omp.h>
#include <set>
#include <vector>

using namespace QutBio;
using namespace std;

// Address of singleton argument table.
Args *arguments;

struct AAClust {
	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using Fasta = FastaSequence;
	using Seq = EncodedFastaSequence;
	using Proto = KmerClusterPrototype;

public:
	static int Run() {
		Params p;
		cerr << arguments->ProgName() << " \\\n" << p << "\n";

		Alphabet * alphabet =  p.alphabet;
		BlosumDifferenceFunction rawDistanceFunction(p.matrix);
		KmerDistanceCache2 distanceFunction(alphabet, &rawDistanceFunction);
		size_t charsPerWord = distanceFunction.CharsPerWord();
		UniformRealRandom rand(p.seed);

		omp_set_num_threads(p.numThreads);

		OMP_TIMER_DECLARE(loadTime);
		OMP_TIMER_DECLARE(clusterTime);

		OMP_TIMER_START(loadTime);

		auto dbSeqs = Load::Fasta(p.fastaFile, p.idIndex, alphabet);
		auto db = Load::Encoded(dbSeqs, -1, alphabet, p.wordLength,
			distanceFunction.CharsPerWord(), alphabet->DefaultSymbol());

		cerr << "AAClust: " << db.size() << " sequences loaded.\n";

		KmerIndex kmerIndex(db, p.wordLength);
		vector<Kmer *> allKmers = kmerIndex.GetKmers();

		OMP_TIMER_END(loadTime);
		OMP_TIMER_START(clusterTime);

		vector<FastaSequence *> protoSeqs;

		auto createPrototype = [&](Kmer *kmer) {
			auto id = String::Format("proto_%6z", protoSeqs.size());
			auto protoSeq = new FastaSequence(id, kmer->Word(), 0, alphabet);
			protoSeqs.push_back(protoSeq);
			auto proto = new KmerClusterPrototype(protoSeq, alphabet, (size_t) p.wordLength, charsPerWord, alphabet->DefaultSymbol());
			return proto;
		};

		map<string, set<string>> homologs;

		if (p.homologFile.size() > 0) {
			homologs = Load::Homologs(p.homologFile);

			cerr << "AAClust: " << homologs.size() << " homolog lists loaded.\n";
		}

		auto assignPurity = [&](KmerClusterPrototype *proto, Cluster *cluster) {
			if (p.homologFile.size() > 0) {
				cluster->AssignPurityMeasure(homologs);
			}
		};

		auto computeDistanceDistribution = [&](KmerClusterPrototype *proto, Cluster *cluster) {
			if (p.computeDistances) {
				cluster->ComputeDistanceDistribution(distanceFunction, p.wordLength, 10);
			}
		};

		auto updateClusterSize = [&](KmerClusterPrototype *proto, Cluster *cluster) {
			proto->Size(proto->Size() + cluster->InstanceCount());
		};

		ofstream protoOut(p.protoOut);
		ofstream clusterOut(p.clusterOut);

		auto saveData = [&](KmerClusterPrototype *proto, Cluster *cluster) {
			if (cluster->kmers.size() > 0) {
				protoOut << proto;
				clusterOut << *cluster;
			}
		};

		//auto showSize = [&](KmerClusterPrototype *proto, Cluster *cluster) {
		//	(cerr << proto->Id() << ": cluster->size() = " << cluster->Size() << "\n").flush();
		//};

		vector<function<void(KmerClusterPrototype *proto, Cluster *cluster)>> process{
			//showSize,
			assignPurity,
			computeDistanceDistribution,
			updateClusterSize,
			saveData
		};

		Cluster::DoExhaustiveIncrementalClustering(
			allKmers,
			p.wordLength,
			p.threshold,
			alphabet->Size(),
			distanceFunction,
			rand,
			p.increment,
			createPrototype,
			process
		);

		OMP_TIMER_END(clusterTime);

#if USE_OMP
		cerr << "Elapsed time loading: " << OMP_TIMER(loadTime) << "\n";
		cerr << "Elapsed time clustering: " << OMP_TIMER(clusterTime) << "\n";
#endif

		return 0;
	}

	struct Params {
		string protoFile;
		string protoOut;
		string fastaFile;
		int numThreads = DEFAULT_THREADS;
		int wordLength;
		int threshold;
		int  seed;
		int idIndex;
		string clusterOut;
		Alphabet * alphabet;
		SimilarityMatrix * matrix;
		string homologFile;
		bool computeDistances = false;
		size_t increment = 1000;

		Params() {
			arguments->Required(protoOut, "protoOut",
				"The name a file that will be overwritten with cluster prototype definitions.");

			arguments->Required(fastaFile, "fastaFile",
				"The name of the FASTA formatted dataset of sequences from which clusters will \n"
				"be derived.");

			arguments->Required(idIndex, "idIndex",
				"The zero-origin location of the sequence Id within the pipe-separated list of \n"
				"metadata fields in the definition line.");

			arguments->Required(seed, "seed", "The random number seed. Integer.");

			arguments->Required(threshold, "threshold", "The threshold distance for cluster inclusion. Integer.");

			arguments->Required(clusterOut, "clusterOut",
				"The name a file that will be overwritten with cluster k-mer definitions.");

			arguments->Required(alphabet, matrix);

			arguments->Optional(numThreads, "numThreads", "The number of OpenMP threads to use in parallel regions.");

			arguments->Optional(computeDistances, "computeDistances",
				"Should we compute an ECD of pairwise distances within clusters. (true|false).");

			arguments->Required(wordLength, "wordLength", "The word length used for kmer tiling.");

			arguments->Optional(protoFile, "protoFile", "The name of a pre-existing list of prototypes, if any.");

			arguments->Optional(homologFile, "homologFile",
				"The name of a homolog file. If present, then a purity score will be \n"
				"calculated and added to the cluster metadata.");

			arguments->Optional(increment, "increment",
				"The number of random kmers to select on each round to become prototypes."
			);

			arguments->Help();

			if (!arguments->Ok()) {
				ostringstream cerr;
				cerr << arguments->ProgName() << ": errors while processing arguments.\nFor help, run: AAClust --help\n";
				throw Exception(cerr.str(), FileAndLine);
			}
		}

		friend ostream & operator<<(ostream & str, const Params & p) {
			str
				<< "--protoFile '" << p.protoFile << "' \\\n"
				<< "--protoOut '" << p.protoOut << "' \\\n"
				<< "--fastaFile '" << p.fastaFile << "' \\\n"
				<< "--numThreads '" << p.numThreads << "' \\\n"
				<< "--wordLength '" << p.wordLength << "' \\\n"
				<< "--threshold '" << p.threshold << "' \\\n"
				<< "--seed '" << p.seed << "' \\\n"
				<< "--idIndex '" << p.idIndex << "' \\\n"
				<< "--clusterOut '" << p.clusterOut << "' \\\n"
				<< "--matrix '" << p.matrix << "' \\\n"
				<< "--homologFile '" << p.homologFile << "' \\\n"
				<< "--increment '" << p.increment << "'\n";

			return str;
		}
	};

};

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

int main(int argc, char *argv[]) {
	try {
		Args args(argc, argv,
			"AAClust: Greedy clustering of Amino Acid kmers by substitution matrix.");

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AAClust::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception &ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}


