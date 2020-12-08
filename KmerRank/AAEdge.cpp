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

#include <bitset>
#include <cstdio>
#include <omp.h>
#include <set>

using namespace QutBio;
using namespace std;

// Address of singleton argument table.
Args *arguments;

struct AAEdge {

	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using Codebook = KmerCodebook<DistanceFunction, Kmer>;
	using Fasta = FastaSequence;
	using Seq = EncodedFastaSequence;
	using Proto = KmerClusterPrototype;

public:
	static int Run() {
		bool ok = true;
		string fastaFile;
		string clusterFile;
		string protoFile;
		string outFile;
		int numThreads = 7;
		int wordLength = 32;
		string merge;
		int idIndex;

		if (arguments->IsDefined("help")) {
			vector<string> text{
				"AAEdge: Gets a list of roughly aligned sub-sequences from a collection of amino acid k-mer clusters.",
				"--help	Gets this text.",
				"--fastaFile	Required. A list of one or more file paths. Each file will be parsed as a FASTA file which contains DNA sequences that have been clustered.",
				"--idIndex	Required. The 0-origin position of the sequence ID field in the pipe-separated definition line.",
				"--clusterFile\t	Required. The name of a file containing the prototypes."
				"--protoFile\t	Required. The name of a file containing the prototypes."
				"--outFile	Required. The name of the output file produced by the program. Each run produces a list of rough pairwise alignments between sequences.",
				"--numThreads	Optional; default value = 7. The number of OpenMP threads to use in parallel regions.",
				"--wordLength	Optional; default value = 32. The word length used for kmer tiling.",
				"--clusterFile	Required. The name of a file that contains a list of k-mer clusters.",
				"--merge	Required. The merge mode used to combine overlapping HSKP (Highly Significant Kmer Pairs). Valid values are:",
				"		consecutive: Merge new edge onto existing edge when both endpoints of the new edge are exact continuations of the previous edge.",
				"		overlapping: Merge new edge onto existing edge when both endpoints of the new edge are within the current extent of the corresponding intervals represented by the existing edge. Some kmers covered by an edge built this way may not be HSKPs.",
				"--matrixId	Optional, default = 62. BLOSUM Matrix ID, one of { 35, 40, 45, 50, 62, 80, 100 }. This is ignored if a custom similarity matrix file is specified.",
				"--matrixFile	Optional. File name for custom similarity matrix. Use this to specify some matrix other than BLOSUM, or if a custom alphabet is in use.",
			};

			for (auto s : text) {
				cerr << s << "\n";
			}

			return 0;
		}

		if (!arguments->Get("fastaFile", fastaFile)) {
			cerr << "Required argument '--fastaFile' not supplied.\n";
			ok = false;
		}

		if (!arguments->Get("clusterFile", clusterFile)) {
			cerr << "Required argument '--clusterFile' not supplied.\n";
			ok = false;
		}

		if (!arguments->Get("protoFile", protoFile)) {
			cerr << "AAEdge: error - required argument '--protoFile' not supplied.\n";
			ok = false;
		}

		if (!arguments->Get("idIndex", idIndex)) {
			cerr << "Required argument '--idIndex' not supplied.\n";
			ok = false;
		}

		if (!arguments->Get("numThreads", numThreads)) {
			cerr << "AAEdge: optional argument '--numThreads' not set; running with default value " << numThreads << ".\n";
		}

		if (!arguments->Get("wordLength", wordLength)) {
			cerr << "AAEdge: optional argument '--wordLength' not set; running with default value " << wordLength << ".\n";
		}

		if (!arguments->Get("merge", merge)) {
			cerr << "AAEdge: required argument '--merge' not provided.\n";
			ok = false;
		}
		else {
			merge = String::ToLowerCase(merge);
			vector<string> validValues{ "consecutive", "overlapping" };
			auto found = std::find(validValues.begin(), validValues.end(), merge);

			if (found == validValues.end()) {
				cerr << "AAEdge: invalid argument '" << merge << "' supplied for '--merge'.\n";
				ok = false;
			}
		}

		if (!arguments->Get("outFile", outFile)) {
			cerr << "AAEdge: required argument '--outFile' not provided.\n";
			ok = false;
		}

		int matrixId;

		if (arguments->IsDefined("matrixId")) {
			if (!(arguments->Get("matrixId", matrixId))) {
				cerr << "Argument 'matrixId' not valid." << endl;
				ok = false;
			}

			vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

			bool found = false;

			for (auto x : matrices) {
				if (x == matrixId) { found = true; }
			}

			if (!found) {
				cerr << "Matrix id not recognised." << endl;
				ok = false;
			}
		}

		if (!ok) {
			cerr << "Invalid command line arguments supplied. For help, run: AAEdge --help\n";
			return 1;
		}

		string matrixFile;
		DistanceType * distanceType = DistanceType::BlosumDistance();

		if (arguments->IsDefined("matrixFile")) {
			arguments->Get("matrixFile", matrixFile);
			distanceType = DistanceType::Custom();
			matrixId = -1;
		}

		bool isCaseSensitive = false;

		if (arguments->IsDefined("isCaseSensitive")) {
			if (!arguments->Get("isCaseSensitive", isCaseSensitive)) {
				cerr << "Invalid data for argument 'isCaseSensitive'." << "\n";
				ok = false;
			}
		}

		SimilarityMatrix * matrix = SimilarityMatrix::GetMatrix(distanceType, matrixId, matrixFile, isCaseSensitive);

		if (ok == false || !matrix) {
			cerr << "Unable to construct similarity matrix. For help, run: AAEdge --help\n";
			return 1;
		}

		Alphabet * alphabet = new Alphabet(matrix);
		BlosumDifferenceFunction rawDistanceFunction(matrix);
		DistanceFunction distanceFunction(alphabet, &rawDistanceFunction);

		omp_set_num_threads(numThreads);

		auto dbSeqs = Load::Fasta(fastaFile, idIndex, alphabet);
		auto db = Load::Encoded(dbSeqs, -1, alphabet, wordLength, distanceFunction.CharsPerWord(), alphabet->DefaultSymbol());
		cerr << arguments->ProgName() << ": " << db.size() << " sequences loaded.\n";

		auto protoSeqs = Load::Fasta(protoFile, 0, alphabet);
		auto protos = Load::Prototypes(protoSeqs, alphabet, wordLength, distanceFunction.CharsPerWord());
		cerr << arguments->ProgName() << ": " << protos.size() << " prototypes loaded.\n";

		Index<Seq> seqIndex(db);
		Index<Proto> protoIndex(protos);
		KmerIndex kmerIndex(db, wordLength);

		auto codebook = Load::Codebook(clusterFile, alphabet, distanceFunction, wordLength, distanceFunction.CharsPerWord(), seqIndex, protoIndex, kmerIndex);

		if (codebook->Size() == 0) {
			ostringstream cerr;
			cerr << "Cluster dataset contains no entries; run terminated.\n";
			throw Exception(cerr.str(), FileAndLine);
		}

		vector<Cluster *> &clusters{ codebook->Codebook() };

		std::sort(clusters.begin(), clusters.end(), [&](Cluster * const& a, Cluster * const& b) {
			return a->kmers.size() > b->kmers.size();
		});

		int i;

		for (i = clusters.size() - 1; i >= 0 && clusters[i]->kmers.size() < 2; i--) {
			delete clusters[i];
			clusters[i] = 0;
		}

		clusters.resize(i + 1);

		vector<Edge> edges;
		vector<Edge> *finalEdges = &edges;

		for (auto cluster : clusters) {
			for (size_t i = 0; i < cluster->kmers.size() - 1; i++) {
				auto & k1 = cluster->kmers[i];
				for (size_t j = i + 1; j < cluster->kmers.size(); j++) {
					auto & k2 = cluster->kmers[j];
					for (auto &instance1 : k1.first->Instances()) {
						for (auto &instance2 : k2.first->Instances()) {
							edges.emplace_back(instance1, instance2, wordLength);
						}
					}
				}
			}
		}

		std::sort(edges.begin(), edges.end());

		vector<Edge> mergedEdges;

		if (merge == "consecutive" || merge == "overlapping") {
			finalEdges = &mergedEdges;
		}

		bool haveCurrent = false;
		Edge current;

		for (auto edge : edges) {
			if (haveCurrent) {
				if (merge == "consecutive") {
					if (edge.id1 == current.id1
						&& edge.id2 == current.id2
						&& edge.start1 == current.stop1 - wordLength + 2
						&& edge.start2 == current.stop2 - wordLength + 2) {
						// Extend by 1.
						current.stop1++;
						current.stop2++;
					}
					else {
						mergedEdges.push_back(current);
						current = edge;
					}
				}
				else if (merge == "overlapping") {
					if (edge.id1 == current.id1
						&& edge.id2 == current.id2
						&& edge.start1 <= current.stop1 + 1
						&& edge.start2 <= current.stop2 + 1) {
						// extend to end of next edge
						current.stop1 = edge.stop1;
						current.stop2 = edge.stop2;
					}
					else {
						mergedEdges.push_back(current);
						current = edge;
					}
				}
				else {
					throw NotImplementedException(FileAndLine);
				}
			}
			else {
				current = edge;
				haveCurrent = true;
			}
		}

		{
			ofstream edgeFile(outFile);
			edgeFile << "Edges\t" << finalEdges->size() << "\n";

			for (auto & edge : *finalEdges) {
				edgeFile << edge << "\n";
			}
		}

		delete codebook;

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
		int retCode = AAEdge::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception &ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}


