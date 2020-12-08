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
#include <cstdio>
#include <stdlib.h>
#include <omp.h>
#include <set>
#include <vector>

using namespace QutBio;
using namespace std;

// Singletons.
Args *arguments;
mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

struct AAClustRecluster {
public:
	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using pCluster = Cluster * ;
	using AssignMode = enum { first = 0, all = 1, nearest = 2 };

	struct Params {
	public:
		string seqFile;
		string protoFile;
		string outFile;
		size_t numThreads = 7;
		size_t wordLength = 0;
		int idIndex = 0;
		bool ok = true;
		DistanceType *distanceType = DistanceType::BlosumDistance();
		Alphabet *alphabet;
		SimilarityMatrix *matrix;
		Distance threshold;
		bool parallelClusters = false;
		AssignMode assignMode;

		Params() {

#define REQ(x,h) arguments->Required(x, #x, h)
#define OPT(x,h) arguments->Optional(x, #x, h)
#define STR(x) STR2(x)
#define STR2(x) #x

			REQ(seqFile,
				"A FASTA file which contains amino acid sequences that will be assigned to clusters \n"
				"centred on the existing set of prototypes.");

			REQ(protoFile, "The name of a file containing the prototypes.");

			REQ(idIndex,
				"The 0-origin position of the sequence ID field in the pipe-separated definition \n"
				"line. Usually this will be 0, but for complicated situations such as my encoding of \n"
				"SwissProt, multiple Id candidates are available.");

			REQ(numThreads, "The number of OpenMP threads to use in parallel regions.");

			arguments->Required(alphabet, matrix);

			REQ(wordLength, "The word length used for kmer tiling.");

			REQ(outFile,
				"The name of a file which will be overwritten with the newly assigned clusters.");

			REQ(threshold, "Positive integer specifying the distance cutoff for assignment of kmers \n"
				"to clusters. A kmer is considered to be a member of the cluster if and only if the \n"
				"distance from kmer to cluster centroid is equal to or less than the threshold \n"
				"distance. The threshold should match that used when the prototypes were initially \n"
				"extracted from the database.");

			arguments->Required(alphabet, matrix);

			OPT(parallelClusters,
				"Do we use the 'parallel k-means' style algorithm, or work cluster-by-cluster. \n"
				"In k-means mode, the clusters are all built concurrently, and nothing \n"
				"(true|false).");

			string assignMode = "all";

			OPT(assignMode,
				"Do we assign kmers to the first cluster, the nearest cluster, or to all clusters \n"
				"within range? (first|nearest|all).");

			String::ToLowerInPlace(assignMode);

			vector<string> assignModes{ "first", "all", "nearest" };
			bool assignModeOk = false;

			for (uint i = 0; i < assignModes.size(); i++) {
				if (assignMode == assignModes[i]) {
					this->assignMode = AssignMode(i);
					assignModeOk = true;
					break;
				}
			}

			if (!assignModeOk) {
				cerr << "Invalid assignMode '" << assignMode << "'\n";
				arguments->Fail();
			}

			if (outFile == seqFile || outFile == protoFile) {
				cerr << arguments->ProgName() << ": Output file " << outFile << " will overwrite one of your input files.\n";
				ok = false;
			}

			if (!ok || !arguments->Ok()) {
				arguments->Help();
				throw Exception("Error processing command line arguments.", FileAndLine);
			}
		}
	};

	using Codebook = KmerCodebook<DistanceFunction, Kmer>;
	using Fasta = FastaSequence;
	using Seq = EncodedFastaSequence;
	using Proto = KmerClusterPrototype;

	static int Run() {
		Params p;

		Alphabet *alphabet = p.alphabet;
		BlosumDifferenceFunction rawDistanceFunction(p.matrix);
		DistanceFunction distanceFunction(alphabet, &rawDistanceFunction);

		omp_set_num_threads(p.numThreads);

		auto dbSeqs = Load::Fasta(p.seqFile, p.idIndex, alphabet);
		auto db = Load::Encoded(dbSeqs, -1, alphabet, p.wordLength, distanceFunction.CharsPerWord(), alphabet->DefaultSymbol());
		cerr << arguments->ProgName() << ": " << db.size() << " sequences loaded.\n";

		auto protoSeqs = Load::Fasta(p.protoFile, 0, alphabet);
		auto protos = Load::Prototypes(protoSeqs, alphabet, p.wordLength, distanceFunction.CharsPerWord());
		cerr << arguments->ProgName() << ": " << protos.size() << " prototypes loaded.\n";

		Index<Seq> seqIndex(db);
		Index<Proto> protoIndex(protos);
		KmerIndex kmerIndex(db, p.wordLength);

		OMP_TIMER_DECLARE(encodeDb);
		OMP_TIMER_START(encodeDb);

		if (p.assignMode == AssignMode::all) {
			Encode(kmerIndex, protos, distanceFunction, p.wordLength, p.threshold, p.outFile);
		}
		else if (p.assignMode == AssignMode::first) {
			cerr << "Please run AAClust to complete this operation.\n";
			return 1;
		}
		else if (p.assignMode == AssignMode::nearest) {
			ClusterNearest(kmerIndex, protos, distanceFunction, p.wordLength, p.threshold, p.outFile);
		}
		else {
			throw Exception("Not implemented", FileAndLine);
		}

		OMP_TIMER_END(encodeDb);

		cerr << "Database encoded in " << OMP_TIMER(encodeDb) << "s.\n";
		Util::Free(protos);
		Util::Free(protoSeqs);
		Util::Free(db);
		Util::Free(dbSeqs);

		// SaveSignatures(db, parms.outFile);
		return 0;
	}

	static void ClusterNearest(
		KmerIndex & kmerIdx,
		vector<Proto *> &protos,
		DistanceFunction &distanceFunction,
		uint K,
		Distance threshold,
		string &outFile
	) {
		const uint C = protos.size();
		vector<Cluster *> clusters;
		vector<Kmer *> kmers = kmerIdx.GetKmers();

		for (auto p : protos) {
			clusters.push_back(new Cluster(p->SingletonKmer(), size_t(0), distanceFunction));
		}

#if USE_OMP
#pragma omp parallel
#endif
		{
			vector<vector<pair<size_t, Distance>>> allocated(C);

			for (auto & a : allocated) {
				a.reserve((kmers.size() + C - 1) / C);
			}

#if USE_OMP
#pragma omp for
#endif
			for (uint i = 0; i < kmers.size(); i++) {
				Distance minDist = numeric_limits<Distance>::max();
				size_t nearestClusterIdx;

				EncodedKmer kmerCode = kmers[i]->PackedEncoding();

				for (uint c = 0; c < C; c++) {
					auto proto = protos[c];
					EncodedKmer centroidCode = proto->SingletonKmer()->PackedEncoding();

					auto dist = distanceFunction(centroidCode, kmerCode, K);

					if (dist <= threshold && dist < minDist) {
						nearestClusterIdx = c;
						minDist = dist;
					}
				}

				if (minDist < numeric_limits<Distance>::max()) {
					allocated[nearestClusterIdx].emplace_back(i, minDist);
				}
			}

#if USE_OMP
#pragma omp critical
#endif
			{
				for (uint c = 0; c < C; c++) {
					for ( auto & i: allocated[c] ) {
						clusters[c]->Add(kmers[i.first], i.second);
					}
				}
			}
		}

		ofstream str(outFile);
		for (auto cluster : clusters) {
			str << cluster;
			delete cluster;
		}
	}

	static void Encode(
		KmerIndex & kmers,
		vector<Proto *> &protos,
		DistanceFunction &distanceFunction,
		uint K,
		Distance threshold,
		string &outFile
	) {
		const uint C = protos.size();

		ofstream str(outFile);

#if USE_OMP
#pragma omp parallel for ordered schedule(static,1)
#endif
		for (uint c = 0; c < C; c++) {
			auto & proto = *protos[c];
			EncodedKmer centroidCode = proto.SingletonKmer()->PackedEncoding();

			Cluster cluster(proto.SingletonKmer(), 0, distanceFunction);

			for (auto & kvp : kmers) {
				auto kmer = kvp.second;
				EncodedKmer kmerCode = kmer->PackedEncoding();

				auto dist = distanceFunction(centroidCode, kmerCode, K);

				if (dist <= threshold) {
					cluster.Add(kmer, dist);
				}
			}

#if USE_OMP
#pragma omp ordered
#endif
			{
				str << cluster;
			}
		}

	}
};

int main(int argc, char *argv[]) {
	try {
		string title = string(argv[0]) + ": Generate a new .clusters file by assigning ( or trying to assign) \n"
			"all k-mers observed in the input sequence to a set of prototypes. K-mer k is assigned \n"
			"to cluster c_i if distance(k,p_i) <= T, where p_i is prototype kmer and T is a real-valued \n"
			"threshold > 0.";

		Args args(argc, argv, title);

		arguments = &args;

		double start_time = omp_get_wtime();
		int retCode = AAClustRecluster::Run();
		double end_time = omp_get_wtime();

		cout << "Elapsed time: " << (end_time - start_time) << "s" << endl;

		return retCode;
	}
	catch (Exception &ex) {
		cerr << ex.File() << "(" << ex.Line() << "): " << ex.what() << "\n";
		return 1;
	}
}
