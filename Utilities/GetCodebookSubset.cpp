#include <map>
#include <string>
#include <vector>

#include <Args.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <KmerCodebook.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <SimilarityMatrix.hpp>
#include <FileUtil.hpp>
#include <DataLoader.hpp>

using namespace std;
using namespace QutBio;

namespace Utilities {
	struct GetCodebookSubset {
		using DistanceFunction = KmerDistanceCache2;
		using Cluster = KmerCluster<DistanceFunction, Kmer>;
		using pCluster = Cluster * ;
		using Codebook = KmerCodebook<DistanceFunction, Kmer>;

		static void Run(int argc, char** argv) {
			Args args(argc, argv);
			Params parms(args);
			Alphabet *alphabet = Alphabets::AA();
			BlosumDifferenceFunction dist(parms.matrix);
			DistanceFunction distanceFunction(alphabet, &dist);

			auto protoSeqs = Load::Fasta(parms.protosIn, 0, alphabet);
			auto protos = Load::Prototypes(protoSeqs, alphabet, parms.kmerLength, distanceFunction.CharsPerWord());

			auto wantedSeqs = Load::Fasta(parms.wantedSeqs, parms.idIndex, alphabet);
			auto wanted = Load::Encoded(wantedSeqs, parms.classIndex, alphabet, parms.kmerLength, distanceFunction.CharsPerWord(), alphabet->DefaultSymbol());
			auto wantedSeqIndex = Load::GetSequenceIndex(wanted);
			auto wantedKmerIndex = Load::GetKmerIndex(wanted, parms.kmerLength);

			vector<pCluster> wantedClusters;
			GetWantedClusters(protos, wantedKmerIndex, distanceFunction, parms.threshold, wantedClusters);

			ofstream pOut(parms.protosOut);
			for (pCluster c : wantedClusters) {
				pOut << &(c->prototype->Sequence());
			}

			ofstream cOut(parms.clustersOut);
			for (pCluster c : wantedClusters) {
				cOut << (*c);
			}
		}

		static void	GetWantedClusters(
			vector<KmerClusterPrototype *> & codebook,
			KmerIndex & wantedKmers,
			DistanceFunction & distance,
			Distance threshold,
			vector<pCluster> & wantedClusters //
		) {
#pragma omp parallel for
			for (uint i = 0; i < codebook.size(); i++) {
				// (cerr << "Processing prototype " << codebook[i]->DefLine() << '\n').flush();
				auto protoKmer = codebook[i]->SingletonKmer();
				auto protoEncoding = protoKmer->PackedEncoding();
				Cluster * newCluster = nullptr;

				for (auto & kmer : wantedKmers) {
					auto encoding = kmer.second->PackedEncoding();
					auto dist = distance(protoEncoding, encoding, protoKmer->Length());

					if (dist <= threshold) {
						if (!newCluster) {
							newCluster = new Cluster(protoKmer, 0, distance);
						}
						newCluster->Add(kmer.second, dist);
					}
				}

				if (newCluster) {
#pragma omp critical
					wantedClusters.push_back(newCluster);
				}
			}
		}

		struct Params {
			bool ok;
			string protosIn, protosOut, clustersOut;
			string wantedSeqs;
			int idIndex, classIndex, kmerLength;
			Distance threshold;
			SimilarityMatrix * matrix;
			Alphabet * alphabet;

			Params(Args & args) {
				ok = true;

				if (!args.Get("protosIn", protosIn)) {
					cerr << "Argument 'protosIn' not supplied.\n";
					ok = false;
				}
				if (!args.Get("protosOut", protosOut)) {
					cerr << "Argument 'protosOut' not supplied.\n";
					ok = false;
				}
				if (!args.Get("clustersOut", clustersOut)) {
					cerr << "Argument 'clustersOut' not supplied.\n";
					ok = false;
				}
				if (!args.Get("wantedSeqs", wantedSeqs)) {
					cerr << "Argument 'wantedSeqs' not supplied.\n";
					ok = false;
				}
				if (!args.Get("idIndex", idIndex)) {
					cerr << "Argument 'idIndex' not supplied.\n";
					ok = false;
				}
				if (!args.Get("classIndex", classIndex)) {
					cerr << "Argument 'classIndex' not supplied.\n";
					ok = false;
				}
				if (!args.Get("kmerLength", kmerLength)) {
					cerr << "Argument 'kmerLength' not supplied.\n";
					ok = false;
				}
				if (!args.Get("threshold", threshold)) {
					cerr << "Argument 'threshold' not supplied.\n";
					ok = false;
				}

				args.Required(alphabet, matrix);

				if (!ok) {
					throw Exception("Invalid arguments.", FileAndLine);
				}
			}

		};

	};
}

using namespace Utilities;

int main(int argc, char** argv) {
	try {
		GetCodebookSubset::Run(argc, argv);
	}
	catch (Exception &ex) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
		return 1;
	}
	catch (runtime_error & err) {
		cerr << "Unhandled exception:" << endl << err.what() << endl;
		return 1;
	}
	return 0;
}

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;
