#pragma once

#include <string>
#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>

#include "QutBio.h"
#include "Random.hpp"
#include "Homologs.hpp"
#include "OpenHausdorff.hpp"
#include "SequenceWrapper.hpp"

using namespace QutBio;

namespace Utilities {

	class GetKmerNNDistributions {
	private:
		struct Parameters {
			string dbFile;
			string homologFile;
			string matrixFile;

			string overallAicFile;
			string overallDistributionFile;

			string posAicFile;
			string posDistributionFile;

			string negAicFile;
			string negDistributionFile;

			uint idIndex = 0;
			uint idDigits = 5;
			int classIndex = -1;
			uint kmerLength = 0;

			DistanceType * dist = DistanceType::HalperinEtAl();
			uint matrixId = 62;

			time_t seed = time(0);

			uint sampleSize = 100;
			uint maxModelSize = 10;

			friend ostream & operator<<(ostream& out, const Parameters & parms) {
				out << "Parameters:" << endl;
				out << "--dbFile " << parms.dbFile << endl;
				out << "--homologFile " << parms.homologFile << endl;
				out << "--matrixFile " << parms.matrixFile << endl;
				out << "--overallAicFile " << parms.overallAicFile << endl;
				out << "--overallDistributionFile " << parms.overallDistributionFile << endl;
				out << "--posAicFile " << parms.posAicFile << endl;
				out << "--posDistributionFile " << parms.posDistributionFile << endl;
				out << "--negAicFile " << parms.negAicFile << endl;
				out << "--negDistributionFile " << parms.negDistributionFile << endl;
				out << "--idIndex " << parms.idIndex << endl;
				out << "--idDigits " << parms.idDigits << endl;
				out << "--classIndex " << parms.classIndex << endl;
				out << "--kmerLength " << parms.kmerLength << endl;
				out << "--dist " << parms.dist->Name() << endl;
				out << "--matrixId " << parms.matrixId << endl;
				out << "--seed " << parms.seed << endl;
				out << "--sampleSize " << parms.sampleSize << endl;
				out << "--maxModelSize " << parms.maxModelSize << endl;
				return out;
			}
		};
	public:

		/**
		Processes a ranking file produced by KmerRank to construct intra-class and
		inter-class empirical distributions of distances.
		The format of the input dataset is (laid out in comma-separated columns):

		query_id,ignored_class_id,subject_id,ignored_class_id,distance,rank,isHomolog

		*/
		static void Run(Args * arguments) {
			Parameters p;
			GetValidatedParameters(*arguments, p);

			cerr << p;

			UniformRealRandom rand((int)p.seed);
			Alphabet * alphabet = Alphabets::AA();
			SimilarityMatrix * matrix = SimilarityMatrix::GetMatrix(alphabet, p.dist, p.matrixId, p.matrixFile);

			auto dbSeqs = Load::Fasta(p.dbFile, p.idIndex, alphabet);
			cerr << arguments->ProgName() << ": " << dbSeqs.size() << " sequences loaded.\n";

			struct Seq : public SequenceWrapper {
				Seq(FastaSequence * base) : SequenceWrapper(base) {}

				vector<Seq *> homologs;

				bool IsHomolog(Seq * other) {
					return find(homologs.begin(), homologs.end(), other) != homologs.end();
				}
			};

			vector<Seq *> db; 
			SequenceWrapper::Wrap(dbSeqs, db);

			Index<Seq> idx(db);

			if (p.homologFile.size() > 0) {
				if (p.homologFile.find_last_of("qrels") == string::npos) {
					ifstream inStream(p.homologFile);
					Homologs::Parse(inStream, idx, idx);
				}
				else {
					FILE * qrelsFile = fopen(p.homologFile.c_str(), "r");

					if (qrelsFile) {
						Homologs::ParseQrels(qrelsFile, idx, idx);
						fclose(qrelsFile);
					}
				}
			}

			vector<Distance> positiveDistances;
			vector<Distance> negativeDistances;

			const size_t N = db.size();

			auto maxKmers = [p](size_t init_, Seq * seq) {
				return std::max<size_t>(init_, seq->KmerCount(p.kmerLength));
			};

			size_t maxKmerCount = std::accumulate(db.begin(), db.end(), 0, maxKmers);
			vector<Distance> rowMinima(maxKmerCount);
			vector<Distance> colMinima(maxKmerCount);

			Projector distCalc(matrix, p.kmerLength, FragmentAggregationMode::HausdorffAverageAverage());

			while (positiveDistances.size() < p.sampleSize || negativeDistances.size() < p.sampleSize) {
				volatile int ix = int(rand() * N);
				volatile int iy = int(rand() * N);

				while (ix == iy) {
					iy = int(rand() * N);
				}

				auto query = db[ix];
				size_t qLen = query->KmerCount(p.kmerLength);

				auto subject = db[iy];
				size_t sLen = subject->KmerCount(p.kmerLength);

				bool homologous = query->IsHomolog(subject);
				distCalc.ComputeDistanceMatrix(
					// *query,
					// *subject,
					query->Sequence().data(),
					subject->Sequence().data(),
					rowMinima,
					colMinima,
					qLen,
					sLen
				);

				auto saveDistance = [&](Distance d) {
					if (homologous) {
						positiveDistances.push_back(d);
					}
					else {
						negativeDistances.push_back(d);
					}
				};

				for (uint i = 0; i < qLen; i++) {
					saveDistance(rowMinima[i]);
				}

				for (uint i = 0; i < sLen; i++) {
					saveDistance(colMinima[i]);
				}
			}

			size_t n = std::min(positiveDistances.size(), negativeDistances.size());

			if (positiveDistances.size() > n) {
				vector<Distance> t;
				Selector s(rand, n, positiveDistances.size());

				for (auto d : positiveDistances) {
					if (s.SelectThis()) t.push_back(d);
				}

				positiveDistances = t;
			}

			if (negativeDistances.size() > n) {
				vector<Distance> t;
				Selector s(rand, n, negativeDistances.size());

				for (auto d : negativeDistances) {
					if (s.SelectThis()) t.push_back(d);
				}

				negativeDistances = t;
			}

			vector<Distance> distances;

			for (size_t i = 0; i < n; i++) {
				distances.push_back(positiveDistances[i]);
				distances.push_back(negativeDistances[i]);
			}

			cerr << "positiveDistances.size " << positiveDistances.size() << endl;
			cerr << "negativeDistances.size " << negativeDistances.size() << endl;
			cerr << "distances.size " << distances.size() << endl;

			ProcessDistribution(p.maxModelSize, p.overallAicFile, p.overallDistributionFile, distances);

			if (positiveDistances.size() > 0) {
				ProcessDistribution(p.maxModelSize, p.posAicFile, p.posDistributionFile, positiveDistances);
				ProcessDistribution(p.maxModelSize, p.negAicFile, p.negDistributionFile, negativeDistances);
			}
		}

		typedef struct {
			GMM1D model;
			double AICc;
		} FittedModel;

		static void EmitDistributions(const string & distributionFileName, const vector<Distance> & distances, vector<FittedModel> & mixtures) {
			ofstream out(distributionFileName);
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

		static void EmitMixtureModels(const string & aicFileName, vector<FittedModel> & mixtures) {
			ofstream out(aicFileName);

			out << "M,AICc" << endl;

			for (uint i = 1; i <= mixtures.size(); i++) {
				FittedModel &m(mixtures[i - 1]);
				out << i << "," << m.AICc << endl;
			}

			out << endl;

			for (uint i = 1; i <= mixtures.size(); i++) {
				out << mixtures[i - 1].model << endl;
			}
		}

		static void ProcessDistribution(size_t maxModelSize, const string & aicFileName, const string & distributionFileName, vector<Distance> & distances) {
			vector<FittedModel> mixtures;
			FitMixtureModels(maxModelSize, distances, mixtures);
			EmitMixtureModels(aicFileName, mixtures);
			sort(distances.begin(), distances.end());
			EmitDistributions(distributionFileName, distances, mixtures);
		}

		/// <summary> Parses and validates the arguments, returning the results in a Parameters object.
		/// </summary>
		/// <param name="arguments"></param>
		/// <returns></returns>

		static void GetValidatedParameters(Args & arguments, Parameters & parms) {
			if (!arguments.Get("dbFile", parms.dbFile)) {
				cerr << "Argument 'dbFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("homologFile", parms.homologFile)) {
				cerr << "Warning: argument 'homologFile' not defined." << endl;
			}

			if (!arguments.Get("overallAicFile", parms.overallAicFile)) {
				cerr << "Argument 'overallAicFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("overallDistributionFile", parms.overallDistributionFile)) {
				cerr << "Argument 'overallDistributionFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("posAicFile", parms.posAicFile)) {
				cerr << "Argument 'posAicFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("posDistributionFile", parms.posDistributionFile)) {
				cerr << "Argument 'posDistributionFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("negAicFile", parms.negAicFile)) {
				cerr << "Argument 'negAicFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("negDistributionFile", parms.negDistributionFile)) {
				cerr << "Argument 'negDistributionFile' not defined." << endl;
				abort();
			}

			if (arguments.IsDefined("idDigits")) {
				arguments.Get("idDigits", parms.idDigits);
			}

			if (arguments.IsDefined("idIndex")) {
				arguments.Get("idIndex", parms.idIndex);
			}

			if (arguments.IsDefined("classIndex")) {
				arguments.Get("classIndex", parms.classIndex);
			};

			if (!(parms.classIndex != (int) parms.idIndex)) {
				cerr << "Argument 'classIndex' must be different from 'idIndex'." << endl;
				abort();
			}

			if (!(arguments.Get("kmerLength", parms.kmerLength))) {
				cerr << "Argument 'kmerLength' not valid." << endl;
				abort();
			}

			if (!(arguments.Get("sampleSize", parms.sampleSize))) {
				cerr << "Argument 'sampleSize' not valid." << endl;
				abort();
			}

			if (!(parms.sampleSize > 0)) {
				cerr << "Argument 'sampleSize' must be greater than zero." << endl;
				abort();
			}

			if (arguments.IsDefined("maxModelSize")) {
				if (!(arguments.Get("maxModelSize", parms.maxModelSize))) {
					cerr << "Argument 'maxModelSize' not valid." << endl;
					abort();
				}
			}

			if (!(arguments.Get("dist", DistanceType::Values(), parms.dist))) {
				cerr << "Argument 'dist' not valid." << endl;
				abort();
			}

			if (arguments.IsDefined("matrixId")) {
				if (!(arguments.Get("matrixId", parms.matrixId))) {
					cerr << "Argument 'matrixId' not valid." << endl;
					abort();
				}
			}
			else {
				parms.matrixId = 62;
			}

			if (parms.dist == DistanceType::Custom()) {
				if (!arguments.IsDefined("matrixFile")) {
					cerr << "Argument 'matrixFile' not supplied for custom similarity matrix." << endl;
					abort();
				}
				else {
					arguments.Get("matrixFile", parms.matrixFile);
				}
			}
			else {
				vector<uint> matrices{ 35, 40, 45, 50, 62, 80, 100 };

				bool found = false;

				for (auto x : matrices) {
					if (x == parms.matrixId) { found = true; }
				}

				if (!found) {
					cerr << "Matrix id not recognised." << endl;
					abort();
				}
			}

			if (arguments.IsDefined("seed")) {
				int seed;

				if (!(arguments.Get("seed", seed))) {
					cerr << "Argument 'seed' not valid." << endl;
					abort();
				}

				if (seed > 0) parms.seed = seed;
			}
		}
	};
}
