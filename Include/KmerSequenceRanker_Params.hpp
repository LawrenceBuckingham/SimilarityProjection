#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include "Args.hpp"
#include "FragmentAggregationMode.hpp"
#include "DistanceType.hpp"
#include "Alphabet.hpp"
#include "AlphabetHelper.hpp"

namespace QutBio {

	struct KmerSequenceRanker_Params {
		string queryFile;
		string dbFile;
		string codebookFile;
		string prototypeFile;
		string matrixFile;

		string rankingFile;
		size_t numThreads = 1;
		int idIndex = 0;
		int classIndex = -1;
		int idDigits = 5;
		int kmerLength = 30;
		int fragLength = 1;
		FragmentAggregationMode * fragMode = FragmentAggregationMode::HausdorffAverage();
		FragmentAggregationMode * kmerMode = FragmentAggregationMode::HausdorffAverage();
		DistanceType * distance = DistanceType::BlosumDistance();
		int matrixId = 62;
		SimilarityMatrix * matrix = 0;
		Alphabet * alphabet = Alphabets::AA();

		Distance thresholdDistance = BAD_DIST;
		Distance defaultDistance = BAD_DIST;

		bool isCaseSensitive = false;

		bool pushKmerDistances = false;

		size_t maxRecords = 1000;

		size_t sampleSize = 0;

		size_t skip = 1;

		int seed = (int)time(NULL);

		string queryIdFile = "";

		/// <summary>
		///		Threshold kmer count for splitting large clusters. 
		///		Leave at -1 to prevent splitting.
		///		Balancing the cluster sizes will not reduce the total 
		///		burden, but it may give better thread utilisation.
		///		A value of 200 seems ok for Test89.
		/// </summary>
		int balancedClusterSize = -1;

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

		KmerSequenceRanker_Params(Args & arguments) {
			bool ok = true;

			if (!arguments.Get("dbFile", dbFile)) {
				cerr << "Argument 'dbFile' not defined." << endl;
				ok = false;
			}

			if (!arguments.Get("queryFile", queryFile)) {
				cerr << "Argument 'queryFile' not defined." << endl;
				ok = false;
			}

			if (arguments.IsDefined("idDigits")) {
				arguments.Get("idDigits", idDigits);
			}

			if (!(idDigits >= 0)) {
				cerr << "Argument 'idDigits' not valid." << endl;
				ok = false;
			}

			if (arguments.IsDefined("idIndex")) {
				arguments.Get("idIndex", idIndex);
			}

			if (!(idIndex >= 0)) {
				cerr << "Argument 'idIndex' not valid." << endl;
				ok = false;
			}

			if (arguments.IsDefined("classIndex")) {
				arguments.Get("classIndex", classIndex);
			};

			if (!(classIndex != idIndex)) {
				cerr << "Argument 'classIndex' must be different from 'idIndex'." << endl;
				ok = false;
			}

			if (!(arguments.Get("fragLength", fragLength) && fragLength > 0)) {
				cerr << "Argument 'fragLength' not valid." << endl;
				ok = false;
			}

			if (!(arguments.Get("fragMode", FragmentAggregationMode::Values(), fragMode) && fragMode != 0)) {
				cerr << "Argument 'fragMode' not valid." << endl;
				ok = false;
			}

			if (!(arguments.Get("kmerMode", FragmentAggregationMode::Values(), kmerMode) && kmerMode != 0)) {
				cerr << "Argument 'kmerMode' not valid." << endl;
				ok = false;
			}

			if (!(arguments.Get("kmerLength", kmerLength))) {
				cerr << "Argument 'kmerLength' not valid." << endl;
				ok = false;
			}

			if (!(kmerLength > 0)) {
				cerr << "Argument 'kmerLength' not valid." << endl;
				ok = false;
			}
			// assert_true( kmerLength <= 9 );

			AlphabetHelper::GetAlphabetAndMatrix(arguments, alphabet, matrix, distance, cerr);

			if (!arguments.Get("rankingFile", rankingFile)) {
				cerr << "Argument 'rankingFile' not defined." << endl;
				ok = false;
			}

			if (!arguments.Get("codebookFile", codebookFile)) {
				if (arguments.IsDefined("prototypeFile")) {
					cerr << "Prototype file is defined but codebook file is not defined.\n";
					ok = false;
				}
				else {
					cerr << "Argument 'codebookFile' not defined.\nCondition treating as warning.\nRelevance must be determined manually." << endl;
				}
			}

			if (!arguments.Get("prototypeFile", prototypeFile)) {
				if (arguments.IsDefined("codebookFile")) {
					cerr << "Codebook file is defined but prototype file is not defined.\n";
					ok = false;
				}
				else {
					cerr << "Argument 'prototypeFile' not defined.\nCondition treating as warning.\nRelevance must be determined manually." << endl;
				}
			}

#define get_opt_arg(v) arguments.GetOptionalArgument( #v, v, ShowHelp )

			get_opt_arg(thresholdDistance);
			get_opt_arg(defaultDistance);
			get_opt_arg(pushKmerDistances);
			get_opt_arg(maxRecords);
			get_opt_arg(sampleSize);
			get_opt_arg(skip);
			get_opt_arg(seed);
			get_opt_arg(queryIdFile);
			get_opt_arg(balancedClusterSize);
			get_opt_arg(numThreads);

#if defined(WANT_DIAGNOSTIC_STREAM) && WANT_DIAGNOSTIC_STREAM
			get_opt_arg(diagnosticsFile);
#endif

			if (skip <= 0) {
				cerr << "Argument 'skip' must be greater than zero." << endl;
			}

			if (seed == -1) {
				seed = (int)time(NULL);
			}

			if (!ok) {
				ShowHelp();
				throw Exception("Invalid arguments.", FileAndLine);
			}
		}

		/// <summary> Displays general help.
		/// </summary>

		static void ShowHelp() {
			string helpText =
				"KmerRank:"  "\n"
				"--------"  "\n"
				"Computes the document ranking of a query dataset under the two - level hierarchical kmer similarity measure."  "\n"
				""  "\n"
				"Required arguments :"  "\n"
				"--------------------"  "\n"
				"--dbFile (fileName) -- The name of a FASTA formatted file containing the database to be queried."  "\n"
				"--queryFile (fileName) -- The name of a FASTA-formatted file containing query sequences. If this"  "\n"
				"    is the same as dbFile, the queries will be sampled uniformly without replacement from the "  "\n"
				"    database." "\n"
				"--kmerLength (uint) -- The number of characters in the kmer tiling."  "\n"
				"--distance (UngappedEdit|HalperinEtAl|Custom) -- The distance measure to use. Custom indicates that a "  "\n"
				"    custom substitution matrix will be supplied."  "\n"
				"--fragLength (uint) -- The fragment length used to partition sequences."  "\n"
				"--fragMode (BestOfBest|Hausdorff|HausdorffAverage) -- The combination mode used to aggregate " "\n"
				"    fragment distances to produce an overall sequence distance."  "\n"
				"--kmerMode (BestOfBest|Hausdorff|HausdorffAverage) -- The combination mode used to aggregate "  "\n"
				"    kmer distances to produce a fragment distance." "\n"
				"--rankingFile (fileName) -- The name of a file which will be overwritten with the CSV-formatted " "\n"
				"    ranking results."  "\n"
				"\n"
				"Optional arguments :"  "\n"
				"--------------------"  "\n"
				"--queryIdFile (fileName) -- Name of a file which, if supplied, will be overwritten with a list of " "\n"
				"    query sequence IDs. This is useful if the queries were sampled from the database, and can be " "\n"
				"    used to extract a suitable subset of the qrels file for accurate Precision-Recall calculations." "\n"
				"--matrixId (35|40|45|50|62|80|100 ) -- The blosum matrix to use, if --distance HalperinEtAl is " "\n"
				"    set. Default value is 62."  "\n"
				"--idIndex (int) -- The zero-origin index position of the ID field in the pipe-separated FASTA " "\n"
				"    definition line. Default value is 0. Each sequence must have a unique Id which will be used"  "\n"
				"    in the output ranking file." "\n"
				"--classIndex (int) -- The zero-origin index position of the class label field in the pipe-separated " "\n"
				"    FASTA definition line. Default value is 1. If no field suits this purpose, use -1."  "\n"
				"    Note that the class label (even if present) is not used in any way by this program."  "\n"
				"--alphabet (AA|DNA|Custom)	-- The alphabet from which symbols are assumed to be drawn. " "\n"
				"    Default is Amino Acid(AA). Use Custom when you want to infer the alphabet from a custom"  "\n"
				"    similarity matrix." "\n"
				"--thresholdDistance (double > 0) -- The threshold distance used to select clusters from the codebook." "\n"
				"--pushKmerDistances (true|false) -- An arcane setting related to some of the exhaustive search algorithms" "\n"
				"    Should I use the \"feed-forward\" mechanism to propagate kmer distances for Hausdorff calculators? " "\n"
				"    Feed-forward is slower if an exhaustive kmer distance table is processed, but it is better able to" "\n"
				"    deal with the situation where the database is scanned kmer-by-kmer rather than sequence by sequence." "\n"
				"--skip (int > 0) -- The sampling rate for the query sequence. Skip = 1 samples every kmer, while " "\n"
				"    skip == 2 samples every second kmer, etc. Default value: 1." "\n"
				"--seed (int) -- Seed for random number generator. " "\n"
				"    Default value (used if seed is missing or negative): time(NULL)." "\n"
				"--maxRecords: (int) -- Optional number of rankings to emit. Default value: 1000." "\n"
				"--sampleSize: (int) -- Optional number of observations to be selected for test. Default value: 1000. " "\n"
				"    Set this to 0 to process the entire query dataset." "\n"
				"--diagnosticsFile: (string) -- Optional name of file which will be overwritten with diagnostic information." "\n"
				;

			cout << helpText;
		}
	};
}
