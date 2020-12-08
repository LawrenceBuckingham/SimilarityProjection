/**
**	Dynamic programming implementation of Similarity projection.
*/

#include "Args.hpp"
#include "KmerSequenceRankerDP.hpp"

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
			Args arguments(argc, argv);

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
			const size_t N = std::min(rankings.size(), maxRecords);

			for (size_t i = 0; i < N; i++) {
				rankingFile << rankings[i] << '\n';
			}
		}

		static void Process(Args & arguments) {
			Params parms{ arguments };

			cout << "--seed " << parms.seed << endl;


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
				exit(1);
			}

			PointerList<FastaSequence> db;

			FastaSequence::ReadSequences(
				db,
				parms.dbFile,
				parms.idIndex,
				parms.classIndex,
				parms.alphabet,
				parms.kmerLength,
				distanceFunction.CharsPerWord(),
				'x',
				FastaSequence::DefaultFactory
			);
			PreprocessDataset(db);

			(cerr << "subject dataset contains " << db.Length() << " sequences." << endl).flush();

			PointerList<FastaSequence> query_;

			if (parms.queryFile != parms.dbFile) {
				FastaSequence::ReadSequences(
					db,
					parms.queryFile,
					parms.idIndex,
					parms.classIndex,
					parms.alphabet,
					parms.kmerLength,
					distanceFunction.CharsPerWord(),
					'x',
					FastaSequence::DefaultFactory
				);
				PreprocessDataset(query_);

				(cerr << "Query dataset contains " << query_.Length() << " sequences." << endl).flush();
			}

			PointerList<FastaSequence> &query = query_.Length() > 0 ? query_ : db;
			vector<FastaSequence *> querySubset;

			if (parms.sampleSize <= 0 || parms.sampleSize >= query.Length()) {
				querySubset = query.Items();
			}
			else {
				UniformRealRandom rand(parms.seed);
				Selector want(rand, parms.sampleSize, query.Length());

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
					fprintf(f, "%s\n", qSeq->Id().c_str());
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
				cerr << "here!\n";
				diagnosticStream = new ofstream(parms.diagnosticsFile);
				ranker.SetDiagnosticStream(diagnosticStream);
			}
#endif
			ranker.RunJob(querySubset, db.Items());
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
		static void PreprocessDataset(PointerList<FastaSequence> & db) {
#pragma omp parallel for
			for (int i = 0; i < (int)db.Length(); i++) {
				auto seq = db[i];
				seq->position = i;
			}
		}

		struct Params {
			string queryFile;
			string dbFile;
			string matrixFile;

			string rankingFile;

			int idIndex = 0;
			int classIndex = -1;
			int idDigits = 5;
			int kmerLength = 30;
			int fragLength = 1;
			FragmentAggregationMode * fragMode = FragmentAggregationMode::HausdorffAverage();
			FragmentAggregationMode * kmerMode = FragmentAggregationMode::HausdorffAverage();
			DistanceType * dist = DistanceType::BlosumDistance();
			int matrixId = 62;
			SimilarityMatrix * matrix;
			Alphabet * alphabet = Alphabet::AA();

			bool isCaseSensitive = false;

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

			Params(Args & arguments) {
				if (!arguments.Get("dbFile", dbFile)) {
					cerr << "Argument 'dbFile' not defined." << endl;
					ShowHelp(); exit(1);
				}

				if (!arguments.Get("queryFile", queryFile)) {
					cerr << "Argument 'queryFile' not defined." << endl;
					ShowHelp(); exit(1);
				}

				if (arguments.IsDefined("idDigits")) {
					arguments.Get("idDigits", idDigits);
				}

				if (!(idDigits >= 0)) {
					cerr << "Argument 'idDigits' not valid." << endl;
					ShowHelp(); exit(1);
				}

				if (arguments.IsDefined("idIndex")) {
					arguments.Get("idIndex", idIndex);
				}

				if (!(idIndex >= 0)) {
					cerr << "Argument 'idIndex' not valid." << endl;
					ShowHelp(); exit(1);
				}

				if (arguments.IsDefined("classIndex")) {
					arguments.Get("classIndex", classIndex);
				};

				if (!(classIndex != idIndex)) {
					cerr << "Argument 'classIndex' must be different from 'idIndex'." << endl;
					ShowHelp(); exit(1);
				}

				if (!(arguments.Get("fragLength", fragLength) && fragLength > 0)) {
					cerr << "Argument 'fragLength' not valid." << endl;
					ShowHelp(); exit(1);
				}

				if (!(arguments.Get("fragMode", FragmentAggregationMode::Values(), fragMode) && fragMode != 0)) {
					cerr << "Argument 'fragMode' not valid." << endl;
					ShowHelp(); exit(1);
				}

				if (!(arguments.Get("kmerMode", FragmentAggregationMode::Values(), kmerMode) && kmerMode != 0)) {
					cerr << "Argument 'kmerMode' not valid." << endl;
					ShowHelp(); exit(1);
				}

				if (!(arguments.Get("kmerLength", kmerLength))) {
					cerr << "Argument 'kmerLength' not valid." << endl;
					ShowHelp(); exit(1);
				}

				if (!(kmerLength > 0)) {
					cerr << "Argument 'kmerLength' not valid." << endl;
					ShowHelp(); exit(1);
				}
				// assert_true( kmerLength <= 9 );

				if (arguments.IsDefined("alphabet")) {
					if (!arguments.Get("alphabet", Alphabet::Values(), alphabet)) {
						cerr << "Unable to parse argument 'alphabet'." << endl;
						ShowHelp(); exit(1);
					}
				}
				else {
					alphabet = Alphabet::AA();
				}

				if (!(arguments.Get("dist", DistanceType::Values(), dist))) {
					cerr << "Argument 'dist' not valid." << endl;
					ShowHelp(); exit(1);
				}

				if (arguments.IsDefined("matrixId")) {
					if (!(arguments.Get("matrixId", matrixId))) {
						cerr << "Argument 'matrixId' not valid." << endl;
						ShowHelp(); exit(1);
					}

					vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

					bool found = false;

					for (auto x : matrices) {
						if (x == matrixId) { found = true; }
					}

					if (!found) {
						cerr << "Matrix id not recognised." << endl;
						ShowHelp(); exit(1);
					}
				}

				bool mustReplaceAlphabet = false;

				if (arguments.IsDefined("matrixFile")) {
					arguments.Get("matrixFile", matrixFile);
					mustReplaceAlphabet = true;
					dist = DistanceType::Custom();
					matrixId = -1;
				}

				arguments.GetOptionalArgument("isCaseSensitive", isCaseSensitive, ShowHelp);

				matrix = SimilarityMatrix::GetMatrix(dist, matrixId, matrixFile, isCaseSensitive);

				if (mustReplaceAlphabet) {
					alphabet = new Alphabet(matrix);
				}

				if (!arguments.Get("rankingFile", rankingFile)) {
					cerr << "Argument 'rankingFile' not defined." << endl;
					ShowHelp(); exit(1);
				}

#define get_opt_arg(v) arguments.GetOptionalArgument( #v, v, ShowHelp )

				get_opt_arg(pushKmerDistances);
				get_opt_arg(maxRecords);
				get_opt_arg(sampleSize);
				get_opt_arg(skip);
				get_opt_arg(seed);
				get_opt_arg(queryIdFile);

#if defined(WANT_DIAGNOSTIC_STREAM) && WANT_DIAGNOSTIC_STREAM
				get_opt_arg(diagnosticsFile);
#endif

				if (skip <= 0) {
					cerr << "Argument 'skip' must be greater than zero." << endl;
				}

				if (seed == -1) {
					seed = (int)time(NULL);
				}
			}

			/// <summary> Displays general help.
			/// </summary>

			static void ShowHelp() {
				string helpText =
					"SimProjDP:"  "\n"
					"----------"  "\n"
					"Computes the document ranking of a query dataset under the two - level hierarchical kmer similarity measure."  "\n"
					""  "\n"
					"Required arguments :"  "\n"
					"--------------------"  "\n"
					"--dbFile (fileName) -- The name of a FASTA formatted file containing the database to be queried."  "\n"
					"--queryFile (fileName) -- The name of a FASTA-formatted file containing query sequences. If this"  "\n"
					"    is the same as dbFile, the queries will be sampled uniformly without replacement from the "  "\n"
					"    database." "\n"
					"--kmerLength (uint) -- The number of characters in the kmer tiling."  "\n"
					"--dist (UngappedEdit|HalperinEtAl|BLOSUM|DES|Custom) -- The distance measure to use. Custom indicates that a "  "\n"
					"    custom substitution matrix will be supplied."  "\n"
					"--fragLength (uint) -- The fragment length used to partition sequences."  "\n"
					"--fragMode (BestOfBest|Hausdorff|HausdorffAverage) -- The combination mode used to aggregate " "\n"
					"    fragment distances to produce an overall sequence distance."  "\n"
					"--kmerMode (BestOfBest|Hausdorff|HausdorffAverage) -- The combination mode used to aggregate "  "\n"
					"    kmer distances to produce a fragment distance." "\n"
					"--rankingFile (fileName) -- The name of a file which will be overwritten with ranking results." "\n"
					"\n"
					"Optional arguments :"  "\n"
					"--------------------"  "\n\n"
					"--queryIdFile (fileName) -- Name of a file which, if supplied, will be overwritten with a list of " "\n"
					"    query sequence IDs. This is useful if the queries were sampled from the database, and can be " "\n"
					"    used to extract a suitable subset of the qrels file for accurate Precision-Recall calculations." "\n\n"

					"--matrixId (35|40|45|50|62|80|100 ) -- The identity of the BLOSUM matrix to use, if --dist is " "\n"
					"    BlosumDistance, HalperinEtAl, or DES. Default value is 62."  "\n\n"

					"--idIndex (int) -- The zero-origin index position of the ID field in the pipe-separated FASTA " "\n"
					"    definition line. Default value is 0. Each sequence must have a unique Id which will be used"  "\n"
					"    in the output ranking file." "\n\n"

					"--classIndex (int) -- The zero-origin index position of the class label field in the pipe-separated " "\n"
					"    FASTA definition line. Default value is 1. If no field suits this purpose, use -1."  "\n"
					"    Note that the class label (even if present) is not used in any way by this program."  "\n\n"

					"--alphabet (AA|DNA|Custom)	-- The alphabet from which symbols are assumed to be drawn. " "\n"
					"    Default is Amino Acid(AA). Use Custom when you want to infer the alphabet from a custom"  "\n"
					"    similarity matrix." "\n\n"

					// TODO: copy this code over to the other version of the program, where it is needed.

					"--thresholdDistance (double > 0) -- The threshold distance used to select clusters from the codebook." "\n"
					"	This setting is applicable only if --dist is " "\n\n"

					"--pushKmerDistances (true|false) -- An arcane setting related to some of the exhaustive search algorithms" "\n"
					"    Should I use the \"feed-forward\" mechanism to propagate kmer distances for Hausdorff calculators? " "\n"
					"    Feed-forward is slower if an exhaustive kmer distance table is processed, but it is better able to" "\n"
					"    deal with the situation where the database is scanned kmer-by-kmer rather than sequence by sequence." "\n\n"

					"--skip (int > 0) -- The sampling rate for the **query** sequence. Skip = 1 samples every kmer, while " "\n"
					"    skip == 2 samples every second kmer, etc. Default value: 1." "\n\n"

					"--seed (int) -- Seed for random number generator. " "\n"
					"    Default value (used if seed is missing or negative): time(NULL)." "\n\n"

					"--maxRecords: (int) -- Optional number of rankings to emit. Default value: 1000." "\n\n"

					"--sampleSize: (int) -- Optional number of observations to be randomly selected for test. Default value: 1000. " "\n"
					"    Set this to 0 to process the entire query dataset." "\n\n"

					"--diagnosticsFile: (string) -- Optional name of file which will be overwritten with diagnostic information." "\n\n"
					;

				cout << helpText;
			}
		};
	};

}

int main(int argc, char** argv) {
	try {
		SimProjDP::Main::Run(argc, argv);
	}
	catch (Exception ex) {
		cerr << "Unhandled exception : " << ex.what() << " - " << ex.File() << "(" << ex.Line() << ")" << endl;
	}
	catch (runtime_error & err) {
		cerr << "Unhandled exception:" << endl << err.what() << endl;
	}
	return 0;
}

mutex FragmentAggregationMode::m;
mutex QutBio::DistanceType::m;

