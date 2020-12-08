#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>
#include <random>

#include "QutBio.h"
#include "Types.hpp"
#include "Homologs.hpp"
#include "DistanceType.hpp"
#include "Random.hpp"
#include "FastaSequence.hpp"
#include "WeibullDistribution.hpp"
#include "NormalDistribution.hpp"

#undef TRON
#include "db.hpp"

using namespace QutBio;
using namespace std;

namespace Utilities
{
typedef unsigned char byte;

class GetKmerTheoreticalDistanceDistributions
{
private:
	struct Parameters
	{
		string dbFile;
		string matrixFile;

		string aicFile;
		string distributionFile;
		string paramFile;

		int idIndex = 0;
		int kmerLength = 0;

		DistanceType *dist = DistanceType::BlosumDistance();
		int matrixId = 62;

		time_t seed = time(0);
		bool isCaseSensitive = false;

		int sampleSize = 10000;
		int maxModelSize = 10;

		pAlphabet alphabet;

		bool ok = true;

		friend ostream &operator<<(ostream &out, const Parameters &parms)
		{
			out << "Parameters:" << endl;
			out << "--dbFile " << parms.dbFile << endl;
			out << "--matrixFile " << parms.matrixFile << endl;
			out << "--aicFile " << parms.aicFile << endl;
			out << "--distributionFile " << parms.distributionFile << endl;
			out << "--idIndex " << parms.idIndex << endl;
			out << "--kmerLength " << parms.kmerLength << endl;
			out << "--dist " << parms.dist->Name() << endl;
			out << "--matrixId " << parms.matrixId << endl;
			out << "--seed " << parms.seed << endl;
			out << "--isCaseSensitive " << parms.isCaseSensitive << endl;
			out << "--sampleSize " << parms.sampleSize << endl;
			out << "--maxModelSize " << parms.maxModelSize << endl;
			out << "--alphabet " << parms.alphabet << endl;
			return out;
		}

		void Help(ostream &out)
		{
			vector<string> help{
				"Parameters:",

				"--dbFile: FileName.",
				"	Path to file containing a FASTA-formatted sequence database.",

				"--aicFile: FileName.",
				"	Path to file that will be overwritten with parameters of fitted Gaussian ",
				"	Mixture Model.",

				"--distributionFile: FileName.",
				"	Path to file that will be overwritten with tabulated empirical and ",
				"	theoretical distribution (CDF) and density (PDF). This file will also ",
				"	include tables containing the exact CDF of the distribution of minimum ",
				"	distances for a range of samples having sizes which are powers of 2, from ",
				"	2 to some maximum value (presently 1024), along with postulated Normal, ",
				"	Weibull, or other analytic distributions which I hope will fit the exact ",
				"	CDF.",
				"",
				"--paramFile: FileName.",
				"	Path to file that will be overwritten with tabulated parameters for the ",
				"	distribution of min{dist(k1_i,k2_i)|i=1..n} for n in {1,2,4,8,16,32,...,",
				"	32768}. Values printed for each n are: mean, standard deviation of exact ",
				"	distribution, and scale and shape parameters of a Weibull approximation.",
				"",
				"--idIndex: int.",
				"	Zero-origin index of pipe-separated field within definition line which ",
				"	contains the sequence id.",

				"--kmerLength: int.",
				"	Length of kmer for distribution calculations.",

				"--dist: BlosumDistance|HalperinEtAl|UngappedEdit|Custom.",
				"	Name of distance function. Best shot is probably BlosumDistance (if you ",
				"	have a similarity matrix) or UngappedEdit (if you just want to count ",
				"	mismatches).",

				"--alphabet: AA|DNA|RNA|DEFAULT|string_to_define_custom_alphabet",
				"	Name of alphabet. AA is the (slightly extended) amino acid alphabet.",
				"	RNA and DNA are the nucleic acid alpahbets, with 'n' representing 'any base'.",
				"	Default is a non-biological alphabet containing the symbols from chr(32)",
				"	up to chr(127). Otherwise, suppply an arbitrary string to define your own alphabet.",

				"--matrixId: int.",
				"	The number of a Blosum matrix. Valid  values are {35, 40, 45, 50, 62, 80, ",
				"	100}. Default = 62.",
				"",
				"--matrixFile: fileName",
				"	Input (optional) custom similarity matrix file. If you don't supply ",
				"	matrixId, then you will need to supply this.",
				"",
				"--isCaseSensitive: true|false.",
				"	Is alphabet case sensitive?",

				"--seed: int.",
				"	Optional seed for random number generator. If not supplied, seed is read ",
				"	from the system clock, timing a resolution of 1 tick per second.",
				"",
				"--sampleSize: int.",
				"	Optional number of points to sample for the empirical CDF (which is ",
				"	intended to act as confirmation of the mixture models). Default ",
				"	value: 10000.",
				"",
				"--maxModelSize: int",
				"	Maximum number of kernels in the Gaussian mixture model.",
			};
			for (auto &s : help)
			{
				cerr << s << "\n\n";
			}
		}

		Parameters(Args &arguments)
		{
			if (arguments.IsDefined("help"))
			{
				Help(cout);
			}

			if (!arguments.Get("dbFile", this->dbFile))
			{
				cerr << "Argument 'dbFile' not defined." << endl;
				ok = false;
			}

			if (arguments.IsDefined("aicFile") && !arguments.Get("aicFile", this->aicFile))
			{
				cerr << "Invalid value for argument 'aicFile'." << endl;
				ok = false;
			}

			if (arguments.IsDefined("distributionFile") && !arguments.Get("distributionFile", this->distributionFile))
			{
				cerr << "Invalid value for argument 'distributionFile'." << endl;
				ok = false;
			}

			if (arguments.IsDefined("paramFile") && !arguments.Get("paramFile", this->paramFile))
			{
				cerr << "Invalid value for argument 'paramFile'." << endl;
				ok = false;
			}

			if (paramFile.length() == 0 && aicFile.length() == 0 && distributionFile.length() == 0)
			{
				cerr << "At least one of 'paramFile', 'aicFile', and 'distributionFile' must be defined." << endl;
				ok = false;
			}

			if (arguments.IsDefined("idIndex"))
			{
				arguments.Get("idIndex", this->idIndex);
			}

			if (!(this->idIndex >= 0))
			{
				cerr << "Argument 'idIndex' not valid." << endl;
				ok = false;
			}

			if (arguments.Get("kmerLength", this->kmerLength))
			{
				if (!(this->kmerLength > 0))
				{
					cerr << "Argument 'kmerLength' not valid." << endl;
					ok = false;
				}
			}
			else
			{
				cerr << "Argument 'kmerLength' not valid." << endl;
				ok = false;
			}

			if (!(arguments.Get("dist", DistanceType::Values(), this->dist)))
			{
				cerr << "Argument 'dist' not valid." << endl;
				ok = false;
			}

			string alpha;

			if (arguments.Get("alphabet", alpha))
			{
				alphabet = Alphabets::ByName(alpha);
			}
			else
			{
				cerr << "Argument 'alphabet' not defined." << endl;
				ok = false;
			}

			if (arguments.IsDefined("matrixId"))
			{
				if (!(arguments.Get("matrixId", this->matrixId)))
				{
					cerr << "Argument 'matrixId' not valid." << endl;
					ok = false;
				}
			}
			else
			{
				this->matrixId = 62;
			}

			if (this->dist == DistanceType::Custom())
			{
				if (!arguments.IsDefined("matrixFile"))
				{
					cerr << "Argument 'matrixFile' not supplied for custom similarity matrix." << endl;
					ok = false;
				}
				else
				{
					arguments.Get("matrixFile", this->matrixFile);
				}
			}
			else
			{
				vector<int> matrices{35, 40, 45, 50, 62, 80, 100};

				bool found = false;

				for (auto x : matrices)
				{
					if (x == this->matrixId)
					{
						found = true;
					}
				}

				if (!found)
				{
					cerr << "Matrix id not recognised." << endl;
					ok = false;
				}
			}

			if (arguments.IsDefined("seed"))
			{
				int seed;

				if (!(arguments.Get("seed", seed)))
				{
					cerr << "Argument 'seed' not valid." << endl;
					ok = false;
				}

				if (seed > 0)
					this->seed = seed;
			}

			if (arguments.IsDefined("isCaseSensitive"))
			{
				if (!(arguments.Get("isCaseSensitive", this->isCaseSensitive)))
				{
					cerr << "Argument 'isCaseSensitive' not valid." << endl;
					ok = false;
				}
			}

			if (arguments.IsDefined("maxModelSize"))
			{
				if (!(arguments.Get("maxModelSize", this->maxModelSize)))
				{
					cerr << "Argument 'maxModelSize' not valid." << endl;
					ok = false;
				}

				if (!(this->maxModelSize > 0))
				{
					cerr << "Argument 'maxModelSize' must be greater than zero." << endl;
					ok = false;
				}
			}

			if (arguments.IsDefined("sampleSize"))
			{
				if (!(arguments.Get("sampleSize", this->sampleSize)))
				{
					cerr << "Argument 'sampleSize' not valid." << endl;
					ok = false;
				}

				if (!(this->sampleSize > 0))
				{
					cerr << "Argument 'sampleSize' must be greater than zero." << endl;
					ok = false;
				}
			}
		}
	};

public:
	struct XPFNE
	{
		double x, p, f, n, e;

		XPFNE(
			double x,
			double p,
			double f,
			double n,
			double e
			//
			) : x(x), p(p), f(f), n(n), e(e)
		{
		}
	};

	/**
		**	Processes a FASTA-formatted sequence file to obtain theoretical and empirical
		**	k-mer distance distributions..
		*/
	static void Run(Args args)
	{
		byte b;
		Parameters parms(args);
		auto alphabet = Alphabets::AA();
		SimilarityMatrix *matrix = SimilarityMatrix::GetMatrix(alphabet, parms.dist, parms.matrixId, parms.matrixFile);

		if (!parms.ok)
		{
			cerr << "Error in command line arguments.\nFor help, GetKmerTheroeticalDistributions --help\n";
			return;
		}

		cerr << parms;

		srand((uint)parms.seed);

		auto db = Load::Fasta(parms.dbFile, parms.idIndex, alphabet);
		cerr << db.size() << " sequences loaded from '" << parms.dbFile << "'" << endl;

		Index<FastaSequence> idx(db);

		auto symbolHistogram = FastaSequence::GetSymbolHistogram(db);
		TRACE;
		function<Distance(Symbol, Symbol)> symbolDistance = [&]( Symbol x, Symbol y) { return matrix->Difference(x, y); };
		TRACE;

		IntegerDistribution rawKmerDist = IntegerDistribution::GetKmerDistanceDistribution<Symbol, Distance>(symbolHistogram, symbolDistance, parms.kmerLength);
		TRACE;

		ofstream *aicFile = parms.aicFile.length() > 0 ? new ofstream(parms.aicFile) : nullptr;
		ofstream *distributionFile = parms.distributionFile.length() > 0 ? new ofstream(parms.distributionFile) : nullptr;
		ofstream *paramFile = parms.paramFile.length() > 0 ? new ofstream(parms.paramFile) : nullptr;
		TRACE;

		if (distributionFile)
			rawKmerDist.Print(*distributionFile, [](double x) { stringstream s; s << setprecision(17) << x; return s.str(); });
		TRACE;

		if (paramFile)
		{
			rawKmerDist.PrintPdf(*paramFile);
		}
		TRACE;

		for (uint n = 2; n <= 1 << 24; n *= 2)
		{
			TRACE;
			auto minDist = rawKmerDist.GetMinimum(n);
			TRACE;
			vector<double> x, F;
			TRACE;
			minDist.TabulateCdf(x, F);
			TRACE;
			WeibullDistribution hope;
			TRACE;
			hope.FitToCdf(x, F);
			TRACE;

			if (paramFile)
			{
				TRACE;
				if (n == 2)
				{
					*paramFile << "Parameters:\n";
					*paramFile << "n\tmu\tsigma\tshape\tscale\n";
				}
				*paramFile << n << "\t" << minDist.Mean() << "\t" << minDist.StdDev() << "\t" << hope.Shape() << "\t" << hope.Scale() << "\n";
			}
			TRACE;

			if (distributionFile)
			{
				*distributionFile << "~~~~~\nn=" << n << "\n~~~~~\n";
				minDist.Print(*distributionFile, [](double x) { stringstream s; s << setprecision(17) << x; return s.str(); });

				*distributionFile << "~~~~~\nn=" << n << "WeibullDistribution(" << hope.Shape() << "," << hope.Scale() << ")\n~~~~~\n";
				*distributionFile << "x\tF\tFitted\n";

				for (size_t i = 0; i < x.size(); i++)
				{
					*distributionFile << x[i] << "\t" << F[i] << "\t" << hope.Cdf(x[i]) << "\n";
				}

				*distributionFile << "~~~~~\nMax absolute error between Sum of m minima of " << n << " distances and its normal approximation\n~~~~~\nm\tmaxError\n";

				IntegerDistribution sumOfM = minDist;

				for (uint m = 2; m <= 1024; m *= 2)
				{
					(cerr << "Processing average of " << m << " minima.\n").flush();

					IntegerDistribution d = sumOfM.AddParallel(sumOfM);
					sumOfM = d;

					ScaledDistribution avg(1.0 / m, d);

					double mu = avg.Mean();
					double sigma = avg.StdDev();
					NormalDistribution norm(mu, sigma);

					double min, max;
					avg.GetSupport(min, max);

					double maxError = 0;

					vector<XPFNE> curves;

					for (double x = min; x <= max; x += (max - min) / 200)
					{
						double p = avg.Pdf(x);
						double F = avg.Cdf(x);
						double N = norm.Cdf(x);
						double error = fabs(F - N);

						if (error > maxError)
						{
							maxError = error;
						}
						curves.emplace_back(x, p, F, N, error);
					}

					(*distributionFile << "Average of minima:\t" << m << "\t"
									   << "Normal approximation:\tmu"
									   << "\t" << mu << "\t"
									   << "sigma"
									   << "\t" << sigma << "\t"
									   << "maxError"
									   << "\t" << maxError << "\n")
						.flush();

					Print(*distributionFile, curves, HighPrecision);
					(*distributionFile << "\n").flush();
				}
			}
			TRACE;

			if (distributionFile)
				(*distributionFile << "\n").flush();
			TRACE;
		}

		if (aicFile && distributionFile)
		{
			// Now, for comparison purposes, do a Monte-Carlo just to make sure we're getting the theoretical
			//	distribution right...
			vector<Distance> distances;
			auto keys = symbolHistogram.GetKeys();
			auto values = symbolHistogram.GetValues();

			IntegerDistribution background(0, (int)(values.size() - 1), values);

			UniformRealRandom random((int32_t)parms.seed);
			Histogram<Symbol> generatedChars;

			for (int i = 0; i < parms.sampleSize; i++)
			{
				Distance d = GetRandomKmerDistance(keys, background, parms.kmerLength, matrix, random, generatedChars);
				distances.push_back(d);
			}

			cerr << "~~~~~~~~~~~~~" << endl;

			symbolHistogram.Print(cerr);

			cerr << "~~~~~~~~~~~~~" << endl;

			generatedChars.Normalise();
			generatedChars.Print(cerr);

			cerr << "~~~~~~~~~~~~~" << endl;

			*distributionFile << endl
							  << endl
							  << endl;

			ProcessDistribution(parms.maxModelSize, *aicFile, *distributionFile, distances);
		}
		TRACE;

		if (aicFile)
			delete aicFile;
		if (distributionFile)
			delete distributionFile;
		if (paramFile)
			delete paramFile;
	}

	static string HighPrecision(double x)
	{
		stringstream s;
		s << setprecision(17) << x;
		return s.str();
	}

	static ostream &Print(
		ostream &out,
		vector<XPFNE> &curves,
		function<string(double)> valFormat)
	{
		out << "x";

		for (auto &c : curves)
		{
			out << "\t" << c.x;
		}

		out << "\nPdf";

		for (auto &c : curves)
		{
			out << "\t" << valFormat(c.p);
		}
		out << "\nCdf";

		for (auto &c : curves)
		{
			out << "\t" << valFormat(c.f);
		}

		out << "\nFitted CDF";

		for (auto &c : curves)
		{
			out << "\t" << valFormat(c.n);
		}

		out << "\nError";

		for (auto &c : curves)
		{
			out << "\t" << fabs(c.e);
		}

		out << "\n";

		return out;
	}

	typedef struct
	{
		GMM1D model;
		double AICc;
	} FittedModel;

	static void EmitEmpiricalDistribution(ostream &out, const vector<Distance> &distances)
	{
		out << "distance";

		for (size_t i = 0; i < distances.size(); i++)
		{
			if (i == 0 || distances[i] != distances[i - 1])
			{
				auto dist = distances[i];
				out << "\t" << dist;
			}
		}

		out << "\nECDF";

		for (size_t i = 0; i < distances.size(); i++)
		{
			if (i == 0 || distances[i] != distances[i - 1])
			{
				out << "\t" << ((double)(i + 1) / distances.size());
			}
		}

		out << "\n";
	}

	static void EmitDistribution(ostream &out, const vector<Distance> &distances, FittedModel &mixture)
	{
		out << "mCDF[" << mixture.model.Size() << "]";

		for (size_t i = 0; i < distances.size(); i++)
		{
			if (i == 0 || distances[i] != distances[i - 1])
			{
				auto dist = distances[i];
				out << "\t" << mixture.model.Cdf(dist);
			}
		}

		out << "\n";
		out << "mPDF[" << mixture.model.Size() << "]";

		for (size_t i = 0; i < distances.size(); i++)
		{
			if (i == 0 || distances[i] != distances[i - 1])
			{
				auto dist = distances[i];
				out << "\t" << mixture.model.Pdf(dist);
			}
		}
		out << "\n";
	}

	static FittedModel FitMixtureModel(size_t modelSize, const vector<Distance> &distances)
	{
		FittedModel m{GMM1D(modelSize), 0};
		m.model.Initialise(distances);
		m.model.Train(distances, 1000, 1e-10, false);
		m.AICc = m.model.AICc(distances);
		cerr << "Gaussian mixture model with " << modelSize << " components: AICc = " << m.AICc << endl;
		return m;
	}

	static void EmitMixtureModel(ofstream &out, FittedModel &mixture)
	{
		out << "Model\nAIC\t" << mixture.AICc << "\n";
		out << mixture.model << "\nEndModel\n";
	}

	static void ProcessDistribution(
		size_t maxModelSize,
		ofstream &aicFile,
		ofstream &distributionFile,
		vector<Distance> &distances)
	{
		vector<Distance> sortedDistances = distances;
		sort(sortedDistances.begin(), sortedDistances.end());
		EmitEmpiricalDistribution(distributionFile, sortedDistances);

		for (uint i = 1; i <= maxModelSize; i++)
		{
			FittedModel mixture = FitMixtureModel(i, distances);
			EmitMixtureModel(aicFile, mixture);
			EmitDistribution(distributionFile, sortedDistances, mixture);
		}
	}

	static Distance GetRandomKmerDistance(
		vector<Symbol> &keys,
		IntegerDistribution &background,
		int kmerLength,
		SimilarityMatrix *matrix,
		UniformRealRandom &random,
		Histogram<Symbol> &generatedChars)
	{
		Distance d = 0;

		for (int i = 0; i < kmerLength; i++)
		{
			double p1 = random();
			int x1 = (int)background.InverseCdf(p1);
			auto c1 = keys[x1];

			double p2 = random();
			int x2 = (int)background.InverseCdf(p2);
			auto c2 = keys[x2];

			d += matrix->Difference(c1, c2);
			generatedChars.Add(c1);
			generatedChars.Add(c2);
		}

		return d;
	}
};
} // namespace AdHoc

mutex QutBio::DistanceType::m;

int main(int argc, char **argv)
{
	Args arguments(argc, argv);

	arguments.Show();

	if (arguments.IsDefined("numThreads"))
	{
		int numThreads;
		if (arguments.Get("numThreads", numThreads))
		{
			omp_set_num_threads(numThreads);
		}
	}

	try
	{
		Utilities::GetKmerTheoreticalDistanceDistributions::Run(arguments);
	}
	catch (Exception &ex)
	{
		string message = ex.what();

		if (message.find("Safe to ignore") != 0)
		{
			cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
		}
	}
	catch (runtime_error &err)
	{
		cerr << "Unhandled exception:" << endl
			 << err.what() << endl;
	}

	return 0;
}
