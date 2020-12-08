#pragma once

#include "QutBio.h"
#include <set>
#include <iomanip>
#include <omp.h>
#include <cstdint>

#if defined(Build_GetKmerDistanceDistribution)
#	include "GetKmerDistanceDistributions.hpp"
#elif defined(Build_GetKmerTheoreticalDistanceDistributions)
#	include "GetKmerTheoreticalDistanceDistributions.hpp"
#else
#include "GetKmerNNDistributions.hpp"
#include "GetKmerNNHomologDistributions.hpp"
#include "PlotMixtureModel.hpp"
#include "GetGmmStats.hpp"
#include "EvalGmm.hpp"
#include "GetPairwiseSimilarityProfiles.hpp"
#include "GetKmerDistanceSample.hpp"
#endif

using namespace QutBio;

namespace Utilities {

	class Main {
	public:
		/// <summary> Main class contains the program entry point: Run().
		/// </summary>
		static void Run(int argc, char ** argv) {
			Args arguments(argc, argv);

			arguments.Show();

			if (arguments.IsDefined("numThreads")) {
				int numThreads;
				if (arguments.Get("numThreads", numThreads)) {
					omp_set_num_threads(numThreads);
				}
			}

			try {
				Process(arguments);
			}
			catch (Exception &e) {
				string message = e.what();

				if (message.find("Safe to ignore") != 0) {
					throw e;
				}
			}
		}

		static void Process(Args & args) {
#if defined(Build_GetKmerDistanceDistribution)
			GetKmerDistanceDistribution::Run(args);
#else
			string program;

			if (!args.Get("program", program)) {
				program == "undefined";
			}

			if (program == "GetIndividualFragDistances") {
				GetIndividualFragmentDistances(args);
			}
			//else if (program == "GetRankingDistributions") {
			//	GetRankingDistributions(args);
			//}
			else if (program == "GetKmerNNDistributions") {
				GetKmerNNDistributions::Run(&args);
			}
			else if (program == "GetKmerNNHomologDistributions") {
				GetKmerNNHomologDistributions::Run(args);
			}
			else if (program == "PlotMixtureModel") {
				PlotMixtureModel::Run(args);
			}
			else if (program == "GetGmmStats") {
				GetGmmStats::Run(args);
			}
			else if (program == "EvalGmm") {
				EvalGmm::Run(args);
			}
			else if (program == "GetPairwiseSimilarityProfiles") {
				GetPairwiseSimilarityProfiles::Run(args);
			}
			else if (program == "GetKmerDistanceSample") {
				GetKmerDistanceSample::Run(args);
			}
			else {
				cerr << "Program [" << program << "] not recognised. Known values are:" << endl;
				cerr << "   GetBitmapRepresentation" << endl;
				cerr << "   GetIndividualFragDistances" << endl;
				cerr << "   GetRankingDistributions" << endl;
				cerr << "   GetKmerDistanceDistribution" << endl;
				cerr << "   GetKmerDistanceSample" << endl;
				cerr << "   GetKmerNNDistributions" << endl;
				cerr << "   GetKmerNNHomologDistributions" << endl;
				cerr << "   GetKmerTheoreticalDistanceDistributions" << endl;
				cerr << "   GetPairwiseSimilarityProfiles" << endl;
				cerr << "   PlotMixtureModel" << endl;
				cerr << "   GetGmmStats" << endl;
				cerr << "   EvalGmm" << endl;
			}
#endif
		}

		/**
			Processes a ranking file produced by KmerRank to construct intra-class and
			inter-class empirical distributions of distances.
			The format of the input dataset is (laid out in comma-separated columns):

			query_id,ignored_class_id,subject_id,ignored_class_id,distance,rank,isHomolog

			*/
#if false
		static void GetRankingDistributions(Args args) {
			string inFile, outFile;

			if (!args.Get("inFile", inFile)) {
				// cerr << "Required argument 'inFile' not found." << endl;
				return;
			}

			if (!args.Get("outFile", outFile)) {
				// cerr << "Required argument 'outFile' not found." << endl;
				return;
			}

			ifstream inStream(inFile);
			CsvReader reader(inStream);

			map<string, vector<double>* > intraClass;
			map<string, vector<double>* > interClass;
			vector<double> globalIntra;
			vector<double> globalInter;

			Func1<Ranking &, bool> process_record = [&](Ranking & record) {
				const string & queryId = record.queryId;

				if (record.query->IsHomolog(record.subject)) {
					auto h = intraClass[queryId];

					if (!h) intraClass[queryId] = h = new vector<double>();

					h->push_back(record.distance);
					globalIntra.push_back(record.distance);
				}
				else {
					auto h = interClass[queryId];

					if (!h) interClass[queryId] = h = new vector<double>();

					h->push_back(record.distance);
					globalInter.push_back(record.distance);
				}


				// return queryId <= "g00050";
				return true;
			};

			auto emit = [](const string & queryId, vector<double> * hist, ostream & out, const char * tag) {
				GMM1D d(10);
				d.Initialise(*hist);
				d.Train(*hist, 1000, 1e-10);

				double lo = d.InverseCdf(1e-5);
				double hi = d.InverseCdf(1 - 1e-5);

				out << "d";

				for (double t = lo; t <= hi; t += (hi - lo) / 100) {
					out << "," << t /* d.Pdf(t) */;
				}

				out << endl;
				out << queryId << "_" << tag;


				for (double t = lo; t <= hi; t += (hi - lo) / 100) {
					out << "," << d.Pdf(t);
				}

				out << endl << endl;
			};

			Action load_complete = [&](void) {
				ofstream outStream(outFile);

				for (auto p : intraClass) {
					const string & queryId = p.first;
					vector<double>* hist = p.second;
					emit(queryId, hist, outStream, "intra");
					emit(queryId, interClass[queryId], outStream, "inter");
				}

				emit("global", &globalIntra, outStream, "intra");
				emit("global", &globalInter, outStream, "inter");
			};

			reader.Stream<Ranking>(process_record, load_complete);
		}
#endif

		static void GetIndividualFragmentDistances(Args args) {
			string inFile, outFile, statsFile;

			if (!args.Get("inFile", inFile)) {
				// cerr << "Required argument 'inFile' not found." << endl;
				return;
			}

			if (!args.Get("outFile", outFile)) {
				// cerr << "Required argument 'outFile' not found." << endl;
				return;
			}

			if (!args.Get("statsFile", statsFile)) {
				// cerr << "Required argument 'statsFile' not found." << endl;
				return;
			}

			ifstream inStream(inFile);
			CsvReader reader(inStream);

			map<string, vector<double>> distance;
			string queryId, qFragId;
			string subjectId, sFragId;

			Func1<vector<string> &, bool> process_record = [&](vector<string> & record) {
				if (record[0] == "query") {
					queryId = record[1];
					(cerr << "queryId: " << queryId << endl).flush();
				}
				else if (record[0] == "subject") {
					subjectId = record[1];
				}
				else {
					string a = record[0];
					qFragId = queryId + "|" + (a.size() == 1 ? "000" : a.size() == 2 ? "00" : a.size() == 3 ? "0" : "") + a;

					string b = record[1];
					sFragId = subjectId + "|" + (b.size() == 1 ? "000" : b.size() == 2 ? "00" : b.size() == 3 ? "0" : "") + b;

					double dist = Double::Parse(record[2]);
					distance[qFragId].push_back(dist);
				}

				// return queryId <= "g00050";
				return true;
			};

			Action load_complete = [&](void) {
				ofstream outStream(outFile);
				ofstream statsStream(statsFile);

				for (auto & p : distance) {
					vector<double> & result = p.second;
					sort(result.begin(), result.end());

					Histogram<double> h;
					h.AddRange(result);
					h.Normalise();
					DiscreteDistribution dist;
					dist.SetPmf(h);

					statsStream << p.first
						<< ",self_match," << dist.Pmf().data.begin()->first
						<< ",mean," << dist.Mean()
						<< ",std," << dist.StdDev()
						<< endl;

					outStream << p.first;

					for (auto d : result) {
						outStream << "," << d;
					}

					outStream << endl;
				}
			};

			reader.StreamRecords(process_record, load_complete);
		}
	};
}


