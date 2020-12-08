#pragma once

#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>

#include "QutBio.h"
#include "Homologs.hpp"
#include "OpenHausdorff.hpp"
#include "SequenceWrapper.hpp"

using namespace QutBio;

namespace Utilities {

	class GetKmerNNHomologDistributions {
	private:
		struct Parameters {
			string dbFile;
			string homologFile;
			string matrixFile;

			string aicFile;
			string distributionFile;

			uint idIndex = 0;
			uint idDigits = 5;
			int classIndex = -1;
			uint kmerLength = 0;

			DistanceType * dist = DistanceType::HalperinEtAl();
			uint matrixId = 62;

			time_t seed = time(0);
			bool isCaseSensitive = false;

			uint sampleSize = 100;
			uint maxModelSize = 10;

			friend ostream & operator<<(ostream& out, const Parameters & parms) {
				out << "Parameters:" << endl;
				out << "--dbFile " << parms.dbFile << endl;
				out << "--homologFile " << parms.homologFile << endl;
				out << "--matrixFile " << parms.matrixFile << endl;
				out << "--aicFile " << parms.aicFile << endl;
				out << "--distributionFile " << parms.distributionFile << endl;
				out << "--idIndex " << parms.idIndex << endl;
				out << "--classIndex " << parms.classIndex << endl;
				out << "--kmerLength " << parms.kmerLength << endl;
				out << "--dist " << parms.dist->Name() << endl;
				out << "--matrixId " << parms.matrixId << endl;
				out << "--seed " << parms.seed << endl;
				out << "--isCaseSensitive " << parms.isCaseSensitive << endl;
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
		static void Run(Args args) {
			Parameters p;
			GetValidatedParameters(args, p);
			cerr << p;
			srand((uint)p.seed);
			auto alphabet = Alphabets::AA();
			SimilarityMatrix * matrix = SimilarityMatrix::GetMatrix(alphabet, p.dist, p.matrixId, p.matrixFile);

			auto dbSeqs = Load::Fasta(p.dbFile, p.idIndex, alphabet);
			cerr << args.ProgName() << ": " << dbSeqs.size() << " sequences loaded.\n";

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

			const size_t N = db.size();

			auto maxKmers = [p](size_t init_, Seq * seq) {
				return std::max<size_t>(init_, seq->KmerCount(p.kmerLength));
			};

			size_t maxKmerCount = std::accumulate(db.begin(), db.end(), 0, maxKmers);

			OpenHausdorff h(
				matrix,
				p.kmerLength,
				FragmentAggregationMode::BestOfBest(),
				FragmentAggregationMode::HausdorffAverageAverage(),
				Alphabets::AA(), // TODO: custom alphabet, derived from similarity matrix.
				32, // TODO: fragLength
				maxKmerCount,
				maxKmerCount
			);

			auto hDist = h.Distances();

			ofstream aicFile(p.aicFile);
			ofstream distributionFile(p.distributionFile);

			for (size_t i = 0; i < p.sampleSize; i++) {
				int ix = rand() % N;
				auto query = db[ix];

				vector<Distance> distances;

				for (auto subject : query->homologs) {
					size_t qLen = query->Sequence().size() - p.kmerLength + 1;
					vector<Distance> qMin(qLen, numeric_limits<Distance>::max());

					size_t sLen = subject->Sequence().size() - p.kmerLength + 1;
					vector<Distance> sMin(sLen, numeric_limits<Distance>::max());

					h.ComputeDistance(*(query->base), *(subject->base));

					for (size_t qPos = 0; qPos < qLen; qPos++) {
						for (size_t sPos = 0; sPos < sLen; sPos++) {
							Distance d = hDist(qPos, sPos);

							if (d < qMin[qPos]) qMin[qPos] = d;
							if (d < sMin[sPos]) sMin[sPos] = d;
						}
					}

					auto saveDistance = [&](Distance d) {
						distances.push_back(d);
					};

					for (uint i = 0; i < qLen; i++) {
						saveDistance(qMin[i]);
					}

					for (uint i = 0; i < sLen; i++) {
						saveDistance(sMin[i]);
					}
				}

				cerr << "distances.size() = " << distances.size() << endl;

				ProcessDistribution(p.maxModelSize, query->IdStr(), aicFile, distributionFile, distances);
			}
		}

		typedef struct {
			GMM1D model;
			double AICc;
		} FittedModel;

		static void EmitDistributions(ostream & out, const string & queryId, const vector<Distance> & distances, vector<FittedModel> & mixtures) {
			out << "#QueryID," << queryId << endl << endl;

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

		static void EmitMixtureModels(ostream & out, const string & queryId, vector<FittedModel> & mixtures) {
			out << "QueryId," << queryId << endl << endl;

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

		static void ProcessDistribution(size_t maxModelSize, const string & queryId, ostream & aicFile, ostream & distributionFile, vector<Distance> & distances) {
			vector<FittedModel> mixtures;
			FitMixtureModels(maxModelSize, distances, mixtures);
			EmitMixtureModels(aicFile, queryId, mixtures);
			sort(distances.begin(), distances.end());
			EmitDistributions(distributionFile, queryId, distances, mixtures);
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
				cerr << "Argument 'homologFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("aicFile", parms.aicFile)) {
				cerr << "Argument 'aicFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("distributionFile", parms.distributionFile)) {
				cerr << "Argument 'distributionFile' not defined." << endl;
				abort();
			}

			if (arguments.IsDefined("idIndex")) {
				arguments.Get("idIndex", parms.idIndex);
			}

			if (arguments.IsDefined("classIndex")) {
				arguments.Get("classIndex", parms.classIndex);
			};

			if (!(arguments.Get("kmerLength", parms.kmerLength))) {
				cerr << "Argument 'kmerLength' not valid." << endl;
				abort();
			}

			if (!(arguments.Get("sampleSize", parms.sampleSize))) {
				cerr << "Argument 'sampleSize' not valid." << endl;
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

			if (arguments.IsDefined("isCaseSensitive")) {
				if (!(arguments.Get("isCaseSensitive", parms.isCaseSensitive))) {
					cerr << "Argument 'isCaseSensitive' not valid." << endl;
					abort();
				}
			}
		}
	};
}
