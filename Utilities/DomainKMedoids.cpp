#include <map>
#include <string>
#include <vector>
#include <omp.h>
#include <unordered_set>

#include <Args.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <KmerCodebook.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <FileUtil.hpp>
#include <KMedoids.hpp>
#include <DataLoader.hpp>

using namespace QutBio;
using namespace std;

namespace Utilities {
	struct DomainKMedoids {
		using DistanceFunction = KmerDistanceCache2;
		using KM = KMedoids<DistanceFunction, Kmer>;
		using Cluster = KM::Cluster;
		using SortMode = KM::SortMode;
		using SelectMode = KM::SelectMode;
		using MedoidMode = KM::MedoidMode;

		static void Run(int argc, char** argv) {
			Args args(argc, argv);
			Params parms(args);
			Alphabet *alphabet = parms.alphabet;
			BlosumDifferenceFunction rawDist(parms.matrix);
			DistanceFunction distance(alphabet, &rawDist);

			omp_set_num_threads(parms.numThreads);

			auto domains = Load::Domains(parms.domains);

			auto dbSeqs = Load::Fasta(parms.db, parms.idIndex, alphabet);
			auto db = Load::Encoded(dbSeqs, parms.classIndex, alphabet, parms.kmerLength, distance.CharsPerWord(), alphabet->DefaultSymbol());
			auto dbIdx = Load::GetSequenceIndex(db);

			vector<const Domain*> domainList;
			for (auto & p : domains) {
				if (parms.wantedDomains.size() == 0 || parms.wantedDomains.find(p.second.pfamId) != parms.wantedDomains.end()) {
					domainList.push_back(&(p.second));
				}
			}

			ofstream protoFile(parms.protos);
			ofstream clusterFile(parms.clusters);
			uint clusterCount = 0;

#pragma omp parallel for schedule(dynamic)
			for (uint i = 0; i < domainList.size(); i++) {
				auto dom = domainList[i];
				vector<Subsequence> domainInstances;
				dom->GetInstances(domainInstances, dbIdx);

				//#pragma omp critical
				//				{
				//					cerr << dom->pfamId << ": domainInstances.size() = " << domainInstances.size() << "\n";
				//				}

				if (domainInstances.size() == 0) continue;

				vector<Kmer *> clusterProtos;
				vector<Cluster *> clusters;

				KM::Partition(domainInstances, clusterProtos, clusters, parms.kmerLength, parms.threshold, parms.seed, *alphabet, distance);

#pragma omp critical
				{
					// cerr << dom->pfamId << ": clusters.size() = " << clusters.size() << "\n";
					for (uint i = 0; i < clusters.size(); i++) {
						Cluster *c = clusters[i];
						Kmer * p = clusterProtos[i];

						// Yuck. KmerClusterPrototype is not thread-safe!!!
						string id = "proto_" + std::to_string(clusterCount++);
						string defline = id + "|" + dom->pfamId + "|size=" + std::to_string(c->Size());
						FastaSequence seq_(defline, p->Word(), 0, alphabet);
						EncodedFastaSequence seq(&seq_, -1, alphabet, parms.kmerLength, distance.CharsPerWord(), alphabet->DefaultSymbol());
						Kmer *newProto = new Kmer(seq, 0, parms.kmerLength);
						c->prototype = newProto;

						(protoFile << &seq).flush();
						(clusterFile << (*c)).flush();

						delete c;
						delete p;
					}
				}
			}
		}

		struct Params {
			bool ok;
			string domains, db, protos, clusters;
			int idIndex, classIndex, kmerLength, seed;
			Alphabet * alphabet;
			SimilarityMatrix *matrix;
			bool isCaseSensitive;
			Distance threshold;
			int numThreads;
			unordered_set<string> wantedDomains;

			Params(Args & args) {
				ok = true;

				if (!args.Get("domains", domains)) {
					cerr << "Argument 'domains' not supplied.\n";
					ok = false;
				}
				if (!args.Get("db", db)) {
					cerr << "Argument 'db' not supplied.\n";
					ok = false;
				}
				if (!args.Get("protos", protos)) {
					cerr << "Argument 'protos' not supplied.\n";
					ok = false;
				}
				if (!args.Get("clusters", clusters)) {
					cerr << "Argument 'clusters' not supplied.\n";
					ok = false;
				}
				if (!args.Get("kmerLength", kmerLength)) {
					cerr << "Argument 'kmerLength' not supplied.\n";
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
				if (!args.Get("isCaseSensitive", isCaseSensitive)) {
					cerr << "Argument 'isCaseSensitive' not supplied.\n";
					ok = false;
				}
				if (!args.Get("threshold", threshold)) {
					cerr << "Argument 'threshold' not supplied.\n";
					ok = false;
				}
				if (!args.Get("seed", seed)) {
					cerr << "Argument 'seed' not supplied.\n";
					ok = false;
				}
				if (!args.Get("numThreads", numThreads)) {
					cerr << "Argument 'numThreads' not supplied.\n";
					ok = false;
				}

				args.Get("wantedDomains", wantedDomains);

				args.Required(alphabet, matrix);

				if (!ok) {
					cerr << "Example:\nDomainKMedoids.exe --wantedDomains PF00001 PF00002 --domains swissprot.domains --db sp500000.faa --protos sp500000.domain.protos --clusters sp500000.domain.clusters --idIndex 2 --classIndex 3 --matrixId 62 --isCaseSensitive false --kmerLength 30 --threshold 305" << '\n';
					throw Exception("Invalid arguments.", FileAndLine);
				}
			}

		};

	};
}

using namespace Utilities;

int main(int argc, char** argv) {
	try {
		DomainKMedoids::Run(argc, argv);
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

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;
