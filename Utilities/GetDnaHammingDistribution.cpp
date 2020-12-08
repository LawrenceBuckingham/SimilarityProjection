#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>
#include <random>

#include "QutBio.h"
#include "Homologs.hpp"
#include "OpenHausdorff.hpp"
#include "Random.hpp"
#include "Exception.hpp"

using namespace QutBio;

class GetDnaHammingDistribution {
private:
	struct Parameters {
		bool ok = true;

		string dbFile;

		string aicFile;
		string distributionFile;

		int idIndex = 0;
		int classIndex = -1;
		int kmerLength = 0;

		time_t seed = time(0);
		bool isCaseSensitive = false;

		int sampleSize = 10000;
		int maxModelSize = 100;

		friend ostream & operator<<(ostream& out, const Parameters & parms) {
			out << "Parameters:\n";
			out << "--dbFile " << parms.dbFile << endl;
			out << "--aicFile " << parms.aicFile << endl;
			out << "--distributionFile " << parms.distributionFile << endl;
			out << "--idIndex " << parms.idIndex << endl;
			out << "--classIndex " << parms.classIndex << endl;
			out << "--kmerLength " << parms.kmerLength << endl;
			out << "--seed " << parms.seed << endl;
			out << "--isCaseSensitive " << parms.isCaseSensitive << endl;
			out << "--sampleSize " << parms.sampleSize << endl;
			out << "--maxModelSize " << parms.maxModelSize << endl;
			return out;
		}

		void Help(ostream& out) {
			out << "GetDnaHammingDistribution Arguments:\n\n";
			
			out << 
				"--help             : " <<
				"Gets this message." <<  "\n\n";

			out << 
				"--dbFile           : " <<
				"The name of a FASTA-formatted DNA sequence (input) file." <<  "\n\n";
			
			out << 
				"--aicFile          : " <<
				"The name of a file that will be overwritten with a list of AIC values\n"
				"                     which correspond to the models written out to the distribution file." <<  "\n\n";
			
			out << 
				"--distributionFile : " << 
				"The name of a file that will be overwritten with the parameters for a\n"
				"                     list of Gaussian mixture models. We can then use the AIC values to\n"
				"                     choose a model for inverse CDF fitting." <<  "\n\n";
			
			out << 
				"--idIndex          : " << 
				"The 0-origin index of the preferred ID field within the pipe-delimited\n"
				"                     FASTA definition lines. Default value = 0." <<  "\n\n";
			
			out << 
				"--classIndex       : " << 
				"The 0-origin index of the preferred class field within the pipe-delimited\n"
				"                     FASTA definition lines. Default value = -1 (corresponds to \"no class info\")." <<  "\n\n";
			
			out << 
				"--kmerLength       : " << 
				"The required kmer length." <<  "\n\n";
			
			out << 
				"--seed             : " << 
				"An integer seed for the random number generator. Default value = time(0)." <<  "\n\n";
			
			out << 
				"--isCaseSensitive  : " << 
				"Boolean indicator whether the alphabet is case sensitive.\n"
				"                     Default value = false." <<  "\n\n";
			
			out << 
				"--sampleSize       : " << 
				"Because I'm paranoid I do a Monte Carlo sample of kmer distances to\n"
				"                     allow me to satisfy myself that the theoretical distribution is OK.\n"
				"                     This parameter is the number of random kmer pairs to use. Default\n"
				"                     value = 10000." <<  "\n\n";
			
			out << 
				"--maxModelSize     : " << 
				"The maximum number of unidimensional Gaussian kernels to add to the\n"
				"                     list of mixture models. Default value = 100." <<  "\n\n";
		}


		/// <summary> Parses and validates the arguments, returning the results in a Parameters object.
		/// </summary>
		/// <param name="arguments"></param>
		/// <returns></returns>

		Parameters(Args & arguments) {
			if (arguments.IsDefined("help")) {
				Help(cerr);
				ok = false;
				return;
			}

			if (!arguments.Get("dbFile", dbFile)) {
				cerr << "Argument 'dbFile' not defined.\n";
				ok = false;
			}

			if (!arguments.Get("aicFile", aicFile)) {
				cerr << "Argument 'aicFile' not defined.\n";
				ok = false;
			}

			if (!arguments.Get("distributionFile", distributionFile)) {
				cerr << "Argument 'distributionFile' not defined.\n";
				ok = false;
			}

			if (arguments.IsDefined("idIndex")) {
				arguments.Get("idIndex", idIndex);
			}

			if (!(idIndex >= 0)) {
				cerr << "Argument 'idIndex' not valid.\n";
				ok = false;
			}

			if (arguments.IsDefined("classIndex")) {
				arguments.Get("classIndex", classIndex);
			};

			if (classIndex == idIndex) {
				cerr << "Argument 'classIndex' must be different from 'idIndex'.\n";
				ok = false;
			}

			if (!(arguments.Get("kmerLength", kmerLength))) {
				cerr << "Argument 'kmerLength' not valid.\n";
				ok = false;
			}

			if (kmerLength <= 0 || kmerLength > 32) {
				cerr << "Argument 'kmerLength' not valid. It should be between 1 and 32 inclusive.\n";
				ok = false;
			}

			if (arguments.IsDefined("seed")) {
				int seed_;

				if (!(arguments.Get("seed", seed_))) {
					cerr << "Argument 'seed' not valid.\n";
					ok = false;
				}

				if (seed_ > 0) seed = seed_;
			}

			if (arguments.IsDefined("isCaseSensitive")) {
				if (!(arguments.Get("isCaseSensitive", isCaseSensitive))) {
					cerr << "Argument 'isCaseSensitive' not valid.\n";
					ok = false;
				}
			}

			if (arguments.IsDefined("maxModelSize")) {
				if (!(arguments.Get("maxModelSize", maxModelSize))) {
					cerr << "Argument 'maxModelSize' not valid.\n";
					ok = false;
				}

				if (!(maxModelSize > 0)) {
					cerr << "Argument 'maxModelSize' must be greater than zero.\n";
					ok = false;
				}
			}

			if (arguments.IsDefined("sampleSize")) {
				if (!(arguments.Get("sampleSize", sampleSize))) {
					cerr << "Argument 'sampleSize' not valid.\n";
					ok = false;
				}

				if (!(sampleSize > 0)) {
					cerr << "Argument 'sampleSize' must be greater than zero.\n";
					ok = false;
				}
			}

			if (!ok) {
				cerr << "For help, use --help\n";
			}
		}

	};
public:

	/**
	Processes a ranking file produced by KmerRank to construct intra-class and
	inter-class empirical distributions of distances.
	The format of the input dataset is (laid out in comma-separated columns):

	query_id,ignored_class_id,subject_id,ignored_class_id,distance,rank,isHomolog

	*/
	static void Run(Args args) {
		Parameters parms(args);
		
		if (!parms.ok) return;

		cerr << parms;

		srand((uint)parms.seed);

		auto db = Load::Fasta(parms.dbFile, parms.idIndex, Alphabets::DNA());
		cerr << db.size() << " sequences loaded from '" << parms.dbFile << "'" << endl;

		Index<FastaSequence> idx(db);
		auto symbolHistogram = FastaSequence::GetSymbolHistogram( db );

		Histogram<Distance> oneMers;
		oneMers.GetOneMerHistogram<Symbol, Symbol>(symbolHistogram, []( Symbol x, Symbol y) { return x == y ? 0 : 1; });

		IntegerDistribution d1(oneMers);
		IntegerDistribution current(d1);

		for (int i = 2; i <= parms.kmerLength; i++) {
			IntegerDistribution d = current.Add(d1);
			current = d;
		}

		ofstream aicFile(parms.aicFile);
		ofstream distributionFile(parms.distributionFile);

		current.Print(distributionFile, [](double x) { stringstream s; s << setprecision(17) << x; return s.str(); });

		// Now, for comparison purposes, do a Monte-Carlo just to make sure we're getting the theoretical 
		//	distribution right...
		vector<Distance> distances;
		vector<Symbol> keys = symbolHistogram.GetKeys();
		vector<double> values = symbolHistogram.GetValues();

		IntegerDistribution background(0, (int)(values.size() - 1), values);
		UniformRealRandom random((int32_t)parms.seed);
		Histogram<Symbol> generatedChars;

		for (int i = 0; i < parms.sampleSize; i++) {
			Distance d = GetRandomKmerDistance(keys, background, parms.kmerLength, random, generatedChars);
			distances.push_back(d);
		}

		cerr << "~~~~~~~~~~~~~\n";

		symbolHistogram.Print(cerr);

		cerr << "~~~~~~~~~~~~~\n";

		generatedChars.Normalise();
		generatedChars.Print(cerr);

		cerr << "~~~~~~~~~~~~~\n";

		distributionFile << endl << endl << endl;

		ProcessDistribution(parms.maxModelSize, aicFile, distributionFile, distances);
	}

	typedef struct {
		GMM1D model;
		double AICc;
	} FittedModel;

	static void EmitDistributions(ostream & out, const vector<Distance> & distances, vector<FittedModel> & mixtures) {
		out << "distance,ECDF";

		for (uint i = 0; i < mixtures.size(); i++) {
			out << ",mCDF[" << mixtures[i].model.Size() << "],mPDF[" << mixtures[i].model.Size() << "]";
		}

		out << endl;

		for (size_t i = 0; i < distances.size(); i++) {
			if (i == 0 || distances[i] != distances[i - 1]) {
				auto dist = distances[i];

				out << dist << "," << ((double)(i + 1) / distances.size());

				for (uint j = 0; j < mixtures.size(); j++) {
					out << "," << mixtures[j].model.Cdf(dist) << "," << mixtures[j].model.Pdf(dist);
				}

				out << endl;
			}
		}

		out << endl << endl << endl;
	}

	static void FitMixtureModels(size_t maxModelSize, const vector<Distance> & distances, vector<FittedModel> & mixtures) {
		for (uint i = 1; i <= maxModelSize; i++) {
			FittedModel m_{ GMM1D(i), 0 };
			mixtures.push_back(m_);
			FittedModel &m(mixtures[i - 1]);
			m.model.Initialise(distances);
			m.model.Train(distances, 1000, 1e-10, false);
			m.AICc = m.model.AICc(distances);
			cerr << "Gaussian mixture model with " << i << " components: AICc = " << m.AICc << endl;
		}
	}

	static void EmitMixtureModels(ostream & out, vector<FittedModel> & mixtures) {
		out << "M,AICc\n";

		for (uint i = 1; i <= mixtures.size(); i++) {
			FittedModel &m(mixtures[i - 1]);
			out << i << "," << m.AICc << endl;
		}

		out << endl;

		for (uint i = 1; i <= mixtures.size(); i++) {
			out << mixtures[i - 1].model << endl;
		}
	}

	static void ProcessDistribution(size_t maxModelSize, ostream & aicFile, ostream & distributionFile, vector<Distance> & distances) {
		vector<FittedModel> mixtures;
		FitMixtureModels(maxModelSize, distances, mixtures);
		EmitMixtureModels(aicFile, mixtures);
		sort(distances.begin(), distances.end());
		EmitDistributions(distributionFile, distances, mixtures);
	}

	static Distance GetRandomKmerDistance(
		vector<Symbol> &keys,
		IntegerDistribution & background,
		int kmerLength,
		UniformRealRandom & random,
		Histogram<Symbol> & generatedChars
	) {
		Distance d = 0;

		for (int i = 0; i < kmerLength; i++) {
			double p1 = random();
			int x1 = (int)background.InverseCdf(p1);
			auto c1 = keys[x1];

			double p2 = random();
			int x2 = (int)background.InverseCdf(p2);
			auto c2 = keys[x2];

			d += ((c1 != c2) || (c1 == 'n') || (c2 == 'n') ? 1 : 0);
			generatedChars.Add(c1);
			generatedChars.Add(c2);
		}

		return d;
	}
};


int main(int argc, char** argv) {
	try {
		Args args(argc, argv);
		GetDnaHammingDistribution::Run(args);
	}
	catch (Exception ex) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")\n";
	}
	catch (runtime_error & err) {
		cerr << "Unhandled exception:\n" << err.what() << endl;
	}
	return 0;
}

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;
