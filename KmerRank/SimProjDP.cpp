/**
**	Dynamic programming implementation of Similarity projection.
*/

#include "Args.hpp"
#include "KmerSequenceRankerDP.hpp"
#include "DataLoader.hpp"

#include <set>
#include <iomanip>
#include <numeric>
#include <cstdio>

using namespace QutBio;
using namespace std;

#undef max

namespace SimProjDP {

	class Main {
	public:
		/// <summary> Main class contains the program entry point: Run().
		/// </summary>
		static void Run(int argc, char ** argv) {
			Args arguments(argc, argv,
				"Uses exact similarity projection to rank sequences from database for each query.");

			arguments.Show();

			if (arguments.IsDefined("numThreads")) {
				int numThreads;
				if (arguments.Get("numThreads", numThreads)) {
					omp_set_num_threads(numThreads);
					cerr << "OMP thread count set to " << numThreads << endl;
				}
			}
			else {
				omp_set_num_threads(1);
			}

			try {
				auto start = omp_get_wtime();
				Process(arguments);
				auto end = omp_get_wtime();
				FILE *f = fopen("KmerRank_time.txt", "a");
				fprintf(f, "Elapsed time: %fs\n", (end - start));
				fclose(f);
			}
			catch (Exception & e) {
				string message = e.what();

				if (message.find("Safe to ignore") == 0) {
					(cerr << e.File() << "(" << e.Line() << "): " << message << "\n").flush();
					exit(1);
				}
			}
		}

	private:

		static void WriteRankingsToFile(
			vector<Ranking> & rankings,
			size_t maxRecords,
			ofstream & rankingFile
		) {
			function<const Ranking *(const Ranking &n)> f = [](const Ranking &n) { return &n; };
			Ranking::SerialiseCompact(rankings, f, rankingFile);
		}

		static void Process(Args & arguments) {
			Params parms{ arguments };
			cerr << arguments.ProgName() << " \\\n" << parms;

			using DistanceFunction = KmerDistanceCache2;
			using KmerType = Kmer;

			BlosumDifferenceFunction dist(parms.matrix);
			DistanceFunction distanceFunction(parms.alphabet, &dist);

			bool filesOk = true;

			if (parms.skip > 1 && parms.fragLength > 1) {
				filesOk = false;
				(cerr << "Present version cannot manage skip > 1 and fragLength > 1.\n").flush();
			}

			if (!File::Exists(parms.dbFile)) {
				filesOk = false;
				(cerr << "Database file " << parms.dbFile << " cannot be opened to read." << endl).flush();
			}

			if (!File::Exists(parms.queryFile)) {
				filesOk = false;
				(cerr << "Database file " << parms.dbFile << " cannot be opened to read." << endl).flush();
			}

			if (!filesOk) {
				throw Exception("Problem found with parameters. See stderr.", FileAndLine);
			}

			using Seq = EncodedFastaSequence;

			auto dbSeqs = Load::Fasta(parms.dbFile, parms.idIndex, parms.alphabet);
			auto db = Load::Encoded(dbSeqs, parms.classIndex, parms.alphabet, parms.kmerLength, distanceFunction.CharsPerWord(), 'x');
			PreprocessDataset(db);
			Index<Seq> seqIdx(db);
			KmerIndex kmerIdx(db, parms.kmerLength);

			(cerr << arguments.ProgName() << ": " << db.size() << " sequences loaded from '" << parms.dbFile << "'.\n").flush();

			vector<FastaSequence *> querySeqs_;
			vector<EncodedFastaSequence *> query_;

			if (parms.queryFile != parms.dbFile) {
				querySeqs_ = Load::Fasta(parms.queryFile, parms.idIndex, parms.alphabet);
				query_ = Load::Encoded(querySeqs_, parms.classIndex, parms.alphabet, parms.kmerLength, distanceFunction.CharsPerWord(), 'x');
				PreprocessDataset(query_);
				(cerr << arguments.ProgName() << ": Query dataset contains " << query_.size() << " sequences." << endl).flush();
			}

			vector<EncodedFastaSequence *> &query = query_.size() > 0 ? query_ : db;
			vector<EncodedFastaSequence *> querySubset;

			if (parms.sampleSize <= 0 || parms.sampleSize >= query.size()) {
				querySubset = query;
			}
			else {
				UniformRealRandom rand(parms.seed);
				Selector want(rand, parms.sampleSize, query.size());

				for (auto seq : query) {
					if (want.SelectThis()) {
						querySubset.push_back(seq);
					}
				}

				assert_equal((size_t)parms.sampleSize, querySubset.size());
			}

			(cerr << "Query subset contains " << querySubset.size() << " sequences." << endl).flush();

			// Write out the query IDs if requested.
			if (parms.queryIdFile.length() > 0) {
				auto f = fopen(parms.queryIdFile.c_str(), "wb");

				for (auto qSeq : querySubset) {
					fprintf(f, "%s\n", qSeq->IdStr().c_str());
				}

				fclose(f);
			}

			ofstream rankingFile(parms.rankingFile);

			Action1<vector<Ranking> &> postProcessRankings = [&](vector<Ranking> & rankings) {
				WriteRankingsToFile(rankings, parms.maxRecords, rankingFile);
			};

			KmerSequenceRankerDP<DistanceFunction, KmerType> ranker(
				parms.matrix,
				parms.kmerMode,
				parms.fragMode,
				parms.alphabet,
				parms.fragLength,
				parms.kmerLength,
				parms.pushKmerDistances,
				distanceFunction,
				parms.skip,
				parms.maxRecords
			);
			ranker.SetQueryComplete(postProcessRankings);

#if defined(WANT_DIAGNOSTIC_STREAM)
			ostream * diagnosticStream = 0;

			if (parms.diagnosticsFile.size() > 0) {
				diagnosticStream = new ofstream(parms.diagnosticsFile);
				ranker.SetDiagnosticStream(diagnosticStream);
			}
#endif
			ranker.RunJob(querySubset, db);
			rankingFile.close();

#if defined(WANT_DIAGNOSTIC_STREAM)
			if (diagnosticStream) {
				diagnosticStream->flush();
				delete diagnosticStream;
			}
#endif
		}

		/**
		*	<summary>
		*		If necessary, pads all sequences out to be as long as the word-size, k.
		*		Encodes the sequences to produce packed word matrices (if we insist that
		*		k be even, then we can collapse the matrices down to a single array, which
		*		would be nice.
		*	</summary>
		*/
		static void PreprocessDataset(vector<EncodedFastaSequence *> & db) {
#pragma omp parallel for
			for (int i = 0; i < (int)db.size(); i++) {
				auto & seq = *db[i];
				seq.position = i;
			}
		}

		struct Params {
			string queryFile;
			string dbFile;
			string matrixFile;

			string rankingFile;

			int idIndex = 0;
			int classIndex = -1;
			uint idDigits = 5;
			uint kmerLength = 30;
			uint fragLength = 1;
			FragmentAggregationMode * fragMode = FragmentAggregationMode::HausdorffAverage();
			FragmentAggregationMode * kmerMode = FragmentAggregationMode::HausdorffAverage();
			DistanceType * dist = DistanceType::BlosumDistance();
			int matrixId = 62;
			SimilarityMatrix * matrix;
			Alphabet * alphabet = Alphabets::AA();

			bool pushKmerDistances = false;

			size_t maxRecords = 1000;

			size_t sampleSize = 1000;

			size_t skip = 1;

			int seed = (int)time(NULL);

			string queryIdFile = "";

#if defined(WANT_DIAGNOSTIC_STREAM) && WANT_DIAGNOSTIC_STREAM
			/// <summary>
			///		Name of a file which will be overwritten with detailed diagnostic data.
			///		This is likely to be a copy of the (row,column)-minima, and possibly the
			///		kmer similarity matrix.
			/// </summary>
			string diagnosticsFile = "";
#endif

			/// <summary> Parses and validates the arguments, returning the results in a Params object.
			/// </summary>
			/// <param name="arguments"></param>
			/// <returns></returns>

			Params(Args & a) {
				a.Required(dbFile, "dbFile",
					"The name of a FASTA formatted file containing the database to be queried.");

				a.Required(queryFile, "queryFile",
					"The name of a FASTA-formatted file containing query sequences. If this"  "\n"
					"is the same as dbFile, the queries will be sampled uniformly without replacement from the "  "\n"
					"database.");

				a.Optional(idDigits, "idDigits", "An old legacy which should be pf no use now.");

				a.Optional(idIndex, "idIndex",
					"The zero-origin index position of the ID field in the pipe-separated FASTA " "\n"
					"definition line. Default value is 0. Each sequence must have a unique Id which will be used"  "\n"
					"in the output ranking file.");

				a.Optional(classIndex, "classIndex",
					"The zero-origin index position of the class label field in the pipe-separated " "\n"
					"FASTA definition line. Default value is 1. If no field suits this purpose, use -1."  "\n"
					"Note that the class label (even if present) is not used in any way by this program.");

				a.Optional(fragLength, "fragLength", "The fragment length used to partition sequences.");

				a.Required("fragMode", FragmentAggregationMode::Values(), fragMode,
					"(BestOfBest|Hausdorff|HausdorffAverage) -- The combination mode used to aggregate " "\n"
					"fragment distances to produce an overall sequence distance.");

				a.Required("kmerMode", FragmentAggregationMode::Values(), kmerMode,
					"(BestOfBest|Hausdorff|HausdorffAverage) -- The combination mode used to aggregate "  "\n"
					"kmer distances to produce a fragment distance.");

				a.Optional(kmerLength, "kmerLength", "The word length for kmer tiling.");

				a.Required("distance", DistanceType::Values(), dist,
					"(UngappedEdit|HalperinEtAl|BLOSUM|DES|Custom) -- The distance measure to use. Custom indicates that a "  "\n"
					"custom substitution matrix will be supplied.");

				a.Required(alphabet, matrix);

				a.Required(rankingFile, "rankingFile",
					"The name of a file which will be overwritten with ranking results.");

#define GET_OPT_ARG(v,h) a.Optional( v, #v, h )

				GET_OPT_ARG(pushKmerDistances, "(true|false) -- An arcane setting related to some of the exhaustive search algorithms" "\n"
					"    Should I use the \"feed-forward\" mechanism to propagate kmer distances for Hausdorff calculators? " "\n"
					"    Feed-forward is slower if an exhaustive kmer distance table is processed, but it is better able to" "\n"
					"    deal with the situation where the database is scanned kmer-by-kmer rather than sequence by sequence.");

				GET_OPT_ARG(maxRecords, "Number of rankings to emit. Default value: 1000");

				GET_OPT_ARG(sampleSize, "number of observations to be randomly selected for test. Default value: 1000. " "\n"
					"    Set this to 0 to process the entire query dataset.");

				GET_OPT_ARG(skip,
					"The sampling rate for the **query** sequence. Skip = 1 samples every kmer, while " "\n"
					"skip == 2 samples every second kmer, etc. Default value: 1.");

				GET_OPT_ARG(seed, "Seed for random number generator. " "\n"
					"Default value (used if seed is missing or negative): time(NULL).");

				GET_OPT_ARG(queryIdFile, "Name of a file which, if supplied, will be overwritten with a list of " "\n"
					"query sequence IDs. This is useful if the queries were sampled from the database, and can be " "\n"
					"used to extract a suitable subset of the qrels file for accurate Precision-Recall calculations.");

#if defined(WANT_DIAGNOSTIC_STREAM) && WANT_DIAGNOSTIC_STREAM
				GET_OPT_ARG(diagnosticsFile, "File in which diagnostic information will be saved.");
#endif

				if (seed == -1) {
					seed = (int)time(NULL);
				}

				if (!a.Ok()) {
					a.Help();
					throw Exception("Error processing arguments.", FileAndLine);
				}
			}

			friend ostream & operator<<(ostream & str, const Params & p) {
				str << "queryFile" << p.queryFile << " \\\n"
					<< "dbFile" << p.dbFile << " \\\n"
					<< "matrixFile" << p.matrixFile << " \\\n"
					<< "rankingFile" << p.rankingFile << " \\\n"
					<< "idIndex" << p.idIndex << " \\\n"
					<< "classIndex" << p.classIndex << " \\\n"
					<< "idDigits" << p.idDigits << " \\\n"
					<< "kmerLength" << p.kmerLength << " \\\n"
					<< "fragLength" << p.fragLength << " \\\n"
					<< "fragMode" << p.fragMode << " \\\n"
					<< "kmerMode" << p.kmerMode << " \\\n"
					<< "distance" << p.dist << " \\\n"
					<< "alphabet" << p.alphabet << " \\\n"
					<< "pushKmerDistances" << p.pushKmerDistances << " \\\n"
					<< "maxRecords" << p.maxRecords << " \\\n"
					<< "sampleSize" << p.sampleSize << " \\\n"
					<< "skip" << p.skip << " \\\n"
					<< "seed" << p.seed << " \\\n"
					<< "queryIdFile" << p.queryIdFile << "\n";

				return str;
			}
		};
	};

}

int main(int argc, char** argv) {
	try {
		SimProjDP::Main::Run(argc, argv);
	}
	catch (Exception &ex) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
	}
	catch (runtime_error & err) {
		cerr << "Unhandled exception:" << endl << err.what() << endl;
	}
	return 0;
}

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

