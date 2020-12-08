/*
 *	Selects a designated number of terms from a vocabulary with MinHash
 *	algorithm.
 */

#include <map>
#include <string>
#include <vector>

#include <Args.hpp>
#include <FastaSequence.hpp>
#include <KmerCodebook.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <SimilarityMatrix.hpp>
#include <FileUtil.hpp>
#include <kNearestNeighbours.hpp>
#include <DataLoader.hpp>

using namespace std;
using namespace QutBio;

struct SelectProtosMinHash {
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

		auto dbSeqs = Load::Fasta(parms.db, parms.idIndex, alphabet);
		auto db = Load::Encoded(dbSeqs, 0, alphabet, parms.kmerLength, distanceFunction.CharsPerWord(), alphabet->DefaultSymbol());
		auto dbIndex = Load::GetSequenceIndex(db);
		auto kmerIndex = Load::GetKmerIndex(db, parms.kmerLength);

		auto protoSeqs = Load::Fasta(parms.protosIn, 0, alphabet);
		auto protos = Load::Prototypes(protoSeqs, alphabet, parms.kmerLength, distanceFunction.CharsPerWord());
		auto protoIndex = Load::GetSequenceIndex(protos);

		auto codebook = Load::Codebook(parms.clustersIn, alphabet, distanceFunction, parms.kmerLength, distanceFunction.CharsPerWord(), dbIndex, protoIndex, kmerIndex);

		KnnVector<Cluster *, size_t> subset(parms.desiredSize, 0);
		Substring::Hash hash;

		for (auto cluster : codebook->codebook) {
			size_t hashCode = hash(cluster->prototype->Substr());

			if (subset.canPush(hashCode)) {
				subset.push(cluster, hashCode);
			}
		}

		ofstream pOut(parms.protosOut);
		ofstream cOut(parms.clustersOut);

		for (pair<size_t, Cluster *> & pair : subset) {
			cOut << *pair.second;
			pOut << &pair.second->prototype->Sequence();
		}
	}

	struct Params {
		bool ok;
		string protosIn, clustersIn, protosOut, clustersOut, db;
		int idIndex, kmerLength, desiredSize;

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
			if (!args.Get("desiredSize", desiredSize)) {
				cerr << "Argument 'desiredSize' not supplied.\n";
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
		SelectProtosMinHash::Run(argc, argv);
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
