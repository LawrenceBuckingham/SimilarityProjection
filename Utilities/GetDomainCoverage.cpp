#include <map>
#include <string>
#include <vector>
#include <cstdio>
#include <thread>
#include <chrono>

#include <Args.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <KmerCodebook.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <DataLoader.hpp>

using namespace std;
using namespace QutBio;

namespace Utilities {
	struct GetDomainCoverage {
		using DistanceFunction = KmerDistanceCache2;
		using Cluster = KmerCluster<DistanceFunction, Kmer>;
		using pCluster = Cluster * ;
		using Codebook = KmerCodebook<DistanceFunction, Kmer>;

		static void Run(int argc, char** argv) {
			Args args(argc, argv);
			Params parms(args);
			Alphabet *alphabet = parms.alphabet;
			BlosumDifferenceFunction dist(parms.matrix);
			DistanceFunction distanceFunction(alphabet, &dist);

			auto allDomains = Load::Domains(parms.domainIn);

			auto dbSeqs = Load::Fasta(parms.seqIn, parms.idIndex, alphabet);
			auto db = Load::Encoded(dbSeqs, parms.classIndex, alphabet, parms.kmerLength, distanceFunction.CharsPerWord(), alphabet->DefaultSymbol());
			auto dbIndex = Load::GetSequenceIndex(db);
			auto kmerIndex = Load::GetKmerIndex(db, parms.kmerLength);

			auto protoSeqs = Load::Fasta(parms.protoIn, 0, alphabet);
			auto prototypes = Load::Prototypes(protoSeqs, alphabet, parms.kmerLength, distanceFunction.CharsPerWord());

			map<string, vector<pair<Domain *, Domain::Entry *>>> entriesPerSequence;
			map<string, int> domSeqCount;

			for (auto & domKVP : allDomains) {
				for (auto & entryKVP : domKVP.second.entries) {
					auto seq = dbIndex.find(entryKVP.first);

					if (seq == dbIndex.end()) continue;

					entriesPerSequence[seq->first].push_back(pair<Domain *, Domain::Entry *>(&domKVP.second, &entryKVP.second));
					domSeqCount[domKVP.first] ++;

					// auto & eee = entriesPerSequence[seq->first];
					// auto & fff = eee.back();
					// fprintf( stderr, "Added entry (%s, %s)\n", fff.first->pfamId.c_str(), fff.second->seqId.c_str() );
				}
			}

			map<string, int> domainProtoCount;
			vector<KmerClusterPrototype*> outOfDomainPrototypes;

			for (auto p : prototypes) {
				//fprintf(stderr, "Processing prototype %s\n", p->DefLine().c_str());
				//fprintf(stderr, "Prototype sequence = [%s]\n", p->Sequence().c_str());

				auto kmer = p->SingletonKmer();
				auto & key = kmer->Substr();
				auto kmerInDatabase = kmerIndex.find(key);

				if (kmerInDatabase == kmerIndex.end()) {
					fprintf(stderr, "Kmer ");
					kmer->Substr().fprint(stderr);
					fprintf(stderr, " not found in database!\n");
				}

				bool found = false;

				for (auto & instance : kmerInDatabase->second->Instances()) {
					auto & seq = instance.sequence;
					auto & seqId = seq.IdStr();

					for (auto & dom_entry : entriesPerSequence[seqId]) {
						for (auto & extent : dom_entry.second->extents) {
							if (extent.begin <= instance.kmerPosition + kmerInDatabase->second->Length() - 1 && instance.kmerPosition <= extent.end) {
								domainProtoCount[dom_entry.first->pfamId] ++;
								found = true;
							}
						}
					}

					if (found) break;
				}

				if (!found) {
					outOfDomainPrototypes.push_back(p);
				}
			}

			FILE * out = fopen(parms.outFile.c_str(), "w");

			if (!out) {
				throw Exception("Error opening output file " + parms.outFile, FileAndLine);
			}

			fprintf(out, "%s\t%s\t%s\n", "PFamId", "Sequence Count", "Prototype count");

			for (auto & domCountKVP : domainProtoCount) {
				fprintf(out, "%s\t%d\t%d\n", domCountKVP.first.c_str(), domSeqCount[domCountKVP.first], domCountKVP.second);
			}

			fclose(out);

			fprintf(stderr, "%zu records written to output file %s\n", domainProtoCount.size(), parms.outFile.c_str());

			FILE * nid = fopen(parms.notInDomain.c_str(), "w");

			if (!nid) {
				throw Exception("Error opening output file " + parms.notInDomain, FileAndLine);
			}

			for (auto p : outOfDomainPrototypes) {
				p->base->fprint(nid);
			}

			fclose(nid);

			fprintf(stderr, "%zu prototypes did not appear to be in a domain. They have been saved in %s\n",
				outOfDomainPrototypes.size(), parms.notInDomain.c_str());
			fflush(stderr);
		}

		struct Params {
			bool ok;
			string domainIn, seqIn, protoIn, outFile, notInDomain;
			int idIndex, classIndex;
			Alphabet * alphabet;
			SimilarityMatrix * matrix;
			int kmerLength;

			Params(Args & args) {
				ok = true;

				if (!args.Get("domainIn", domainIn)) {
					cerr << "Argument 'domainIn' not supplied.\n";
					ok = false;
				}
				if (!args.Get("protoIn", protoIn)) {
					cerr << "Argument 'protoIn' not supplied.\n";
					ok = false;
				}
				if (!args.Get("seqIn", seqIn)) {
					cerr << "Argument 'seqIn' not supplied.\n";
					ok = false;
				}
				if (!args.Get("outFile", outFile)) {
					cerr << "Argument 'outFile' not supplied.\n";
					ok = false;
				}
				if (!args.Get("notInDomain", notInDomain)) {
					cerr << "Argument 'notInDomain' not supplied.\n";
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
		GetDomainCoverage::Run(argc, argv);
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
