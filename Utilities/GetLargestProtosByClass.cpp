/*
 *	Extracts a subset of a set of clusters based on class labels.
 *	The program reads prototypes and clusters and selects the K
 *	largest clusters (based on prototype size property) from each
 *	class and saves them as a new file.
 *	Prototypes are expected to have id in field 0 and class label
 *	in field 1 of the FASTA defline.
 */

#include <map>
#include <string>
#include <vector>
#include <cstdio>

#include <Args.hpp>
#include <Alphabet.hpp>
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

struct GetLargestProtosByClass {
	using DistanceFunction = KmerDistanceCache2;
	using Cluster = KmerCluster<DistanceFunction, Kmer>;
	using pCluster = Cluster * ;
	using Codebook = KmerCodebook<DistanceFunction, Kmer>;

	static void Run(int argc, char** argv) {
		Args args(argc, argv);
		Params parms(args);

		//	The alphabet and distance function are dummies required to satisfy
		//	the rather top-heavy API for reading a codebook. That will be fixed 
		//	one day.
		Alphabet *alphabet = Alphabets::AA();
		BlosumDifferenceFunction dist(SimilarityMatrix::Blosum62());
		DistanceFunction distanceFunction(alphabet, &dist);

		auto dbSeqs = FastaSequence::Read(parms.db, parms.idIndex, alphabet);
		auto db = Load::Encoded(dbSeqs, parms.classIndex, alphabet, parms.kmerLength, 
			distanceFunction.CharsPerWord(), alphabet->DefaultSymbol());

		Index<EncodedFastaSequence> dbIndex(db);
		KmerIndex kmerIndex(db, parms.kmerLength);

		auto protoSeqs = Load::Fasta(parms.protosIn, 0, alphabet);
		auto protos = Load::Prototypes( protoSeqs, alphabet, parms.kmerLength, distanceFunction.CharsPerWord());
		Index<KmerClusterPrototype> protoIndex(protos);

		FILE * f = fopen(parms.clustersIn.c_str(), "r");

		if (!f) {
			ostringstream cerr;
			cerr << "Error reading codebook from '" << parms.clustersIn << "'\n";
			throw Exception(cerr.str(), FileAndLine);
		}

		Codebook * codebook = new Codebook(
			alphabet, distanceFunction, distanceFunction.CharsPerWord(), parms.kmerLength, dbIndex,
			protoIndex, kmerIndex, f
		);

		fclose(f);

		map<string, vector<pKmerClusterPrototype>> seqFamilies;

		for (auto p : protos) {
			seqFamilies[p->ClassLabel()].push_back(p);
		}

		map<string, pCluster> clusters;

		for (auto c : codebook->codebook) {
			clusters[c->prototype->Sequence().IdStr()] = c;
		}

		ofstream pOut(parms.protosOut);
		ofstream cOut(parms.clustersOut);

		for (auto & pair : seqFamilies) {
			auto & protosPerFamily = pair.second;

			std::sort(protosPerFamily.begin(), protosPerFamily.end(),
				[](const pKmerClusterPrototype & lhs, const pKmerClusterPrototype & rhs) {
				return lhs->Size() > rhs->Size();
			}
			);

			for (uint i = 0; i < parms.protosPerClass && i < protosPerFamily.size(); i++) {
				auto p = protosPerFamily[i];
				pOut << p;
				cOut << *clusters[p->IdStr()];
			}
		}
	}

	struct Params {
		bool ok;
		string protosIn, clustersIn, protosOut, clustersOut, db;
		uint idIndex, protosPerClass, kmerLength;
		int classIndex;

		Params(Args & args) {
			ok = true;

			if (!args.Get("db", db)) {
				cerr << "Argument 'db' not supplied.\n";
				ok = false;
			}
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
			if (!args.Get("clustersIn", clustersIn)) {
				cerr << "Argument 'clustersIn' not supplied.\n";
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
			if (!args.Get("protosPerClass", protosPerClass)) {
				cerr << "Argument 'protosPerClass' not supplied.\n";
				ok = false;
			}
			if (!args.Get("kmerLength", kmerLength)) {
				cerr << "Argument 'kmerLength' not supplied.\n";
				ok = false;
			}

			if (!ok) {
				throw Exception("Invalid arguments.", FileAndLine);
			}
		}

	};

};

int main(int argc, char** argv) {
	try {
		GetLargestProtosByClass::Run(argc, argv);
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
