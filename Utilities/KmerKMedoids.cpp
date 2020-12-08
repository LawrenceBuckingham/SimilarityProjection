#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include <Args.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <KmerIndex.hpp>
#include <Random.hpp>
#include <KmerDistanceCache.hpp>
#include <KmerCluster.hpp>
#include <FileUtil.hpp>
#include <KMedoids.hpp>
#include <DataLoader.hpp>

// file://a:\phd\Notes\KmerKMedoids.html

using namespace std;
using namespace QutBio;

namespace Utilities {
	struct KmerKMedoids {
		using DistanceFunction = KmerDistanceCache2;
		using KM = KMedoids<DistanceFunction, Kmer>;
		using Cluster = KM::Cluster;
		using SortMode = KM::SortMode;
		using SelectMode = KM::SelectMode;
		using MedoidMode = KM::MedoidMode;

		struct Params {
			bool ok;
			string seqIn, protoOut, clusterOut;
			int idIndex, kmerLength, iterations, seed, trials;
			SortMode sortMode;
			SelectMode selectMode;
			Distance threshold;
			Alphabet * alphabet;
			SimilarityMatrix * matrix;
			MedoidMode medoidMode;
			uint minMedditSize = 1000;

			Params(Args & args) {
				ok = true;

				if (!args.Get("seqIn", seqIn)) {
					cerr << "Argument 'seqIn' not supplied.\n";
					ok = false;
				}

				if (!args.Get("protoOut", protoOut)) {
					cerr << "Argument 'protoOut' not supplied.\n";
					ok = false;
				}

				if (!args.Get("clusterOut", clusterOut)) {
					cerr << "Argument 'clusterOut' not supplied.\n";
					ok = false;
				}

				if (!args.Get("idIndex", idIndex)) {
					cerr << "Argument 'idIndex' not supplied.\n";
					ok = false;
				}

				if (!args.Get("iterations", iterations)) {
					cerr << "Argument 'iterations' not supplied.\n";
					ok = false;
				}

				if (!args.Get("kmerLength", kmerLength)) {
					cerr << "Argument 'kmerLength' not supplied.\n";
					ok = false;
				}

				if (!args.Get("seed", seed)) {
					cerr << "Argument 'seed' not supplied.\n";
					ok = false;
				}

				if (!args.Get("threshold", threshold)) {
					cerr << "Argument 'threshold' not supplied.\n";
					ok = false;
				}

				if (!args.Get("trials", trials)) {
					cerr << "Argument 'trials' not supplied.\n";
					ok = false;
				}

				int t = 0;

				if (!args.Get("sortMode", t)) {
					cerr << "Argument 'sortMode' not supplied.\n";
					ok = false;
				}

				if (t < 1 || t > 3) {
					cerr << "Invalid sortMode: value values are SortRandom = 1, SortLongestFirst = 2, SortShortestFirst = 3.\n";
					ok = false;
				}
				else {
					sortMode = (SortMode)t;
				}

				t = 0;

				if (!args.Get("selectMode", t)) {
					cerr << "Argument 'selectMode' not supplied.\n";
					ok = false;
				}

				if (t < 1 || t > 3) {
					cerr << "Invalid selectMode: value values are SelectGreedy = 1, SelectNearest = 2.\n";
					ok = false;
				}
				else {
					selectMode = (SelectMode)t;
				}

				t = 0;

				if (!args.Get("medoidMode", t)) {
					cerr << "Argument 'medoidMode' not supplied.\n";
					ok = false;
				}

				if (t < (int) KM::MedoidMode::MedoidMin || t > (int) KM::MedoidMode::MedoidMax) {
					cerr << "Invalid medoidMode: value values are MedoidBruteForce = 1, MedoidMeddit = 2.\n";
					ok = false;
				}
				else {
					medoidMode = (MedoidMode)t;
				}

				if (args.IsDefined("minMedditSize") && args.Get("minMedditSize", minMedditSize)) {
					if (minMedditSize < 100) {
						cerr << "Warning: minMedditSize must be greater than or equal to 100.\n";
						minMedditSize = 100;
					}
				}

				args.Required(alphabet, matrix);

				if (!ok) {
					throw Exception("Invalid arguments.", FileAndLine);
				}
			}
		};

		static void Run(int argc, char** argv) {
			Args args(argc, argv);
			Params parms(args);
			Alphabet * alphabet = parms.alphabet;

			BlosumDifferenceFunction rawDist(parms.matrix);
			DistanceFunction distance(alphabet, &rawDist);

			auto dbSeqs = Load::Fasta(parms.seqIn, parms.idIndex, alphabet);
			auto db = Load::Encoded(dbSeqs, 0, alphabet, parms.kmerLength, distance.CharsPerWord(), alphabet->DefaultSymbol());

			UniformIntRandom<uint> rand(parms.seed, 0, db.size() - 1);

			vector<Kmer *> clusterProtos;
			vector<Cluster *> clusters;

			vector<Subsequence> initialSubseqs(db.size());

			auto toSubsequence = [](const pEncodedFastaSequence seq) {
				return Subsequence{ seq, 0, seq->Length() };
			};

			std::transform(db.begin(), db.end(), initialSubseqs.begin(), toSubsequence);

			KM::Partition(
				initialSubseqs, 
				clusterProtos, 
				clusters, 
				parms.kmerLength, 
				parms.threshold, 
				parms.seed, 
				*alphabet, 
				distance, 
				parms.trials, 
				parms.iterations, 
				parms.sortMode, 
				parms.selectMode, 
				parms.medoidMode, 
				parms.minMedditSize
			);

			ofstream cOut(parms.clusterOut);
			ofstream pOut(parms.protoOut);

			for (uint i = 0; i < clusters.size(); i++) {
				Cluster *c = clusters[i];
				Kmer * p = clusterProtos[i];

				// Yuck. KmerClusterPrototype is not thread-safe!!!
				string id = "proto_" + std::to_string(i);
				string defline = id + "|size=" + std::to_string(c->Size());
				FastaSequence seq_(defline, p->Word(), 0, alphabet);
				EncodedFastaSequence seq(&seq_, -1, alphabet, parms.kmerLength, distance.CharsPerWord(), alphabet->DefaultSymbol());
				Kmer * newProto = new Kmer(seq, 0, parms.kmerLength);
				c->prototype = newProto;

				(pOut << &seq).flush();
				(cOut << (*c)).flush();

				delete c;
				delete p;
			}
		}

	};
}

using namespace Utilities;

std::mutex DistanceType::m;

int main(int argc, char** argv) {
	try {
		KmerKMedoids::Run(argc, argv);
	}
	catch (Exception & ex) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
		return 1;
	}
	catch (runtime_error & err) {
		cerr << "Unhandled exception:" << endl << err.what() << endl;
		return 1;
	}
	return 0;
}

