#include <map>
#include <string>
#include <vector>

#include <Args.hpp>
#include <Domain.hpp>
#include <FastaSequence.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <DataLoader.hpp>

// a:\Phd\Notes\ExtractDomainSubsequences.html

using namespace std;
using namespace QutBio;

namespace Utilities {
	struct ExtractDomainSubsequences {
		static void Run(int argc, char** argv) {
			Args args(argc, argv);
			Params parms(args);
			Alphabet *alphabet = Alphabets::AA();

			auto dbSeqs = Load::Fasta(parms.seqIn, parms.idIndex, alphabet);
			auto db = Load::Encoded(dbSeqs, -1, alphabet, parms.kmerLength, 1, alphabet->DefaultSymbol());
			auto seqIndex = Load::GetSequenceIndex(db);

			map<string, Domain> domains;
			ifstream domainFile(parms.domainIn);
			Domain::Load(domainFile, domains);

			ofstream out(parms.seqOut);

			for (auto & domainId : parms.domainIds) {
				Domain & dom = domains[domainId];

				for (auto & entryRec : dom.entries) {
					const string & seqId = entryRec.first;

					auto seqPos = seqIndex.find(seqId);

					if (seqPos == seqIndex.end()) continue;

					uint extentCtr = 0;
					uint numExtents = entryRec.second.extents.size();

					for (auto & extent : entryRec.second.extents) {
						SubVector<Symbol> subseq( & (seqPos->second->Sequence()), extent.begin, extent.end - extent.begin + 1);
						out << ">" << domainId << "_" << seqId;

						if (numExtents > 1) {
							out << "_" << (extentCtr++);
						}

						out << "|len=" << (extent.end - extent.begin + 1) << "\n" << alphabet->Decode(subseq) << "\n";
					}
				}
			}
		}

		struct Params {
			bool ok;
			string seqIn, domainIn, seqOut;
			int idIndex, kmerLength;
			vector<string> domainIds;

			Params(Args & args) {
				ok = true;

				if (!args.Get("seqIn", seqIn)) {
					cerr << "Argument 'seqIn' not supplied.\n";
					ok = false;
				}

				if (!args.Get("domainIn", domainIn)) {
					cerr << "Argument 'domainIn' not supplied.\n";
					ok = false;
				}

				if (!args.Get("seqOut", seqOut)) {
					cerr << "Argument 'seqOut' not supplied.\n";
					ok = false;
				}

				if (!args.Get("idIndex", idIndex)) {
					cerr << "Argument 'idIndex' not supplied.\n";
					ok = false;
				}

				if (!args.Get("kmerLength", kmerLength)) {
					cerr << "Argument 'kmerLength' not supplied.\n";
					ok = false;
				}

				if (!args.Get("domainIds", domainIds)) {
					cerr << "Argument 'domainIds' not supplied.\n";
					ok = false;
				}

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
		ExtractDomainSubsequences::Run(argc, argv);
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

