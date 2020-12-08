
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
#include "FileUtil.hpp"

#include <bitset>
#include <cstdio>
#include <omp.h>
#include <set>

using namespace QutBio;
using namespace std;

// Address of singleton argument table.
Args *arguments;

struct AAClusterFirst {

	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using pCluster = Cluster * ;
	using Codebook = KmerCodebook<DistanceFunction, Kmer>;
	using Fasta = FastaSequence;
	using Seq = EncodedFastaSequence;
	using Proto = KmerClusterPrototype;

	struct Params {
		bool ok = true;

		string fastaFile;
		string clusterFile;
		string protoFile;
		string clusterOut;
		string protoOut;
		size_t numThreads = 7;
		size_t wordLength = 32;
		int idIndex = 0;
		size_t numClusters = 0;
		DistanceType * distanceType = DistanceType::BlosumDistance();
		Alphabet * alphabet;
		SimilarityMatrix * matrix;

		Params() {
			if (arguments->IsDefined("help")) {
				vector<string> text{
					"AAClusterFirst: Gets the ${numClusters} largest clusters from a kmer codebook.",
					"--help	Gets this text.",
					"--fastaFile    Required. A list of one or more file paths. Each file will be parsed as a FASTA file which contains DNA sequences that have been clustered.",
					"--clusterFile    Required. The name of a file that contains a list of k-mer clusters.",
					"--protoFile      Required. The name of a file containing the prototypes."
					"--numClusters  Required. The number of clusters to select fromt he codebook.",
					"--clusterOut   Required. The name of the cluster output file.",
					"--protoOut     Required. The name of the prototype output file.",
					"--idIndex      Required. The 0-origin position of the sequence ID field in the pipe-separated definition line.",
					"--numThreads   Optional; default value = 7. The number of OpenMP threads to use in parallel regions.",
					"--wordLength   Optional; default value = 32. The word length used for kmer tiling.",
					"--matrixId     Optional, default = 62. BLOSUM Matrix ID, one of { 35, 40, 45, 50, 62, 80, 100 }. This is ignored if a custom similarity matrix file is specified.",
					"--matrixFile   Optional. File name for custom similarity matrix. Use this to specify some matrix other than BLOSUM, or if a custom alphabet is in use.",
				};

				for (auto s : text) {
					cerr << s << "\n";
				}
			}

			if (!arguments->Get("fastaFile", fastaFile)) {
				cerr << arguments->ProgName() << ": error - required argument '--fastaFile' not supplied.\n";
				ok = false;
			}

			if (!arguments->Get("clusterFile", clusterFile)) {
				cerr << arguments->ProgName() << ": error - required argument '--clusterFile' not supplied.\n";
				ok = false;
			}

			if (!arguments->Get("protoFile", protoFile)) {
				cerr << arguments->ProgName() << ": error - required argument '--protoFile' not supplied.\n";
				ok = false;
			}

			if (!arguments->Get("idIndex", idIndex)) {
				cerr << arguments->ProgName() << ": error - required argument '--idIndex' not supplied.\n";
				ok = false;
			}

			if (!arguments->Get("numClusters", numClusters)) {
				cerr << arguments->ProgName() << ": error - required argument '--numClusters' not supplied.\n";
				ok = false;
			}

			if (!arguments->Get("numThreads", numThreads)) {
				cerr << arguments->ProgName() << ": note - optional argument '--numThreads' not set; running with default value " << numThreads << ".\n";
			}

			if (!arguments->Get("wordLength", wordLength)) {
				cerr << arguments->ProgName() << ": note - optional argument '--wordLength' not set; running with default value " << wordLength << ".\n";
			}

			if (!arguments->Get("clusterOut", clusterOut)) {
				cerr << arguments->ProgName() << ": note - required argument '--clusterOut' not provided.\n";
				ok = false;
			}

			if (!arguments->Get("protoOut", protoOut)) {
				cerr << arguments->ProgName() << ": note - required argument '--protoOut' not provided.\n";
				ok = false;
			}

			if (!ok) {
				ostringstream cerr;
				cerr << "Invalid command line arguments supplied. For help, run: " << arguments->ProgName() << " --help\n";
				throw Exception(cerr.str(), FileAndLine);
			}

			arguments->Required(alphabet,matrix);

			if (clusterOut == clusterFile || clusterOut == fastaFile || clusterOut == protoFile) {
				ostringstream cerr;
				cerr << "AASample: Output file " << clusterOut << " will overwrite one of your input files.\n";
				throw Exception(cerr.str(), FileAndLine);
			}

			if (protoOut == clusterFile || protoOut == fastaFile || protoOut == protoFile) {
				ostringstream cerr;
				cerr << "AASample: Output file " << protoOut << " will overwrite one of your input files.\n";
				throw Exception(cerr.str(), FileAndLine);
			}
		}
	};

public:
	static int Run() {
		Params p;

		if (!p.ok) {
			throw Exception("Invalid arguments. For help, ask for '--help'.", FileAndLine);
		}

		Alphabet * alphabet = p.alphabet;
		BlosumDifferenceFunction rawDistanceFunction(p.matrix);
		DistanceFunction distanceFunction(alphabet, &rawDistanceFunction);

		omp_set_num_threads(p.numThreads);

		auto dbSeqs = Load::Fasta(p.fastaFile, p.idIndex, alphabet);
		auto db = Load::Encoded(dbSeqs, -1, alphabet, p.wordLength, distanceFunction.CharsPerWord(), 'x');
		cerr << arguments->ProgName() << ": " << db.size() << " sequences loaded.\n";

		auto protoSeqs = Load::Fasta(p.protoFile, 0, alphabet);
		auto protos = Load::Prototypes(protoSeqs, alphabet, p.wordLength, distanceFunction.CharsPerWord());
		cerr << arguments->ProgName() << ": " << protos.size() << " prototypes loaded.\n";

		Index<Seq> seqIndex(db);
		KmerIndex kmerIndex(db, p.wordLength);
		Index<Proto> protoIndex(protos);


		FILE * f = fopen(p.clusterFile.c_str(), "r");

		if (!f) {
			ostringstream cerr;
			(cerr << arguments->ProgName() << ": " << p.clusterFile << " cannot be opened for reading.\n").flush();
			throw Exception(cerr.str(), FileAndLine);
		}

		Codebook * codebook = new Codebook(
			alphabet, distanceFunction, distanceFunction.CharsPerWord(),
			p.wordLength, seqIndex, protoIndex, kmerIndex, f);

		fclose(f);

		if (codebook->Size() == 0) {
			ostringstream cerr;
			cerr << "Cluster dataset contains no entries; run terminated.\n";
			throw Exception(cerr.str(), FileAndLine);
		}

		vector<pCluster> &clusters{ codebook->Codebook() };

		auto descendingClusterSize = [](const pCluster & lhs, const pCluster & rhs) {
			return lhs->InstanceCount() > rhs->InstanceCount();
		};

		std::sort(clusters.begin(), clusters.end(), descendingClusterSize);

		(cerr << "selecting largest " << p.numClusters << " clusters from " << clusters.size() << "\n").flush();

		vector<pCluster> clusterSubset;
		vector<Proto *> protoSubset;

		for (auto cluster : clusters) {
			if (clusterSubset.size() >= p.numClusters) break;

			clusterSubset.push_back(cluster);
			protoSubset.push_back((Proto *)& cluster->prototype->Sequence());
		}

		ofstream cOut(p.clusterOut);
		for (auto c : clusterSubset) cOut << (*c);

		ofstream pOut(p.protoOut);
		for (auto p : protoSubset) pOut << p;

		Util::Free(protos);
		Util::Free(protoSeqs);
		Util::Free(db);
		Util::Free(dbSeqs);

		return 0;
	}
};

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

int main(int argc, char *argv[]) {
	try {
		Args args(argc, argv);

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AAClusterFirst::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception &ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}


