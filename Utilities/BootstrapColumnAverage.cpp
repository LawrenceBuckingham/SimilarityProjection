/**
 *	bootstrap_column_average: gets the average value of
 *	each column in a tab-delimited table, along with bootstrap
 *	confidence intervals.
 */

#if __cplusplus < 201103L
#undef _cplusplus
#define __cplusplus 201703L
#endif

#include <limits>
#include <vector>
#include <fstream>

#include "CsvIO.hpp"
#include "Args.hpp"
#include "Random.hpp"

using namespace std;
using namespace QutBio;

namespace Bootstrap {
	struct Params {
		uint firstColumn = 0;
		uint lastColumn = numeric_limits<int>::max();
		uint bootstrap = 200;
		double lowerLimit = 0.05;
		double upperLimit = 0.95;
		int skipLeadingRows = 1;
		int skipTrailingRows = 1;
		string inFile, outFile;

		static void Help() {
			vector<string> help{
				"Compute the expected value, plus bootstrap confidence interval,",
				"of a range of numeric columns in a tab-separated tabular text file."
				"",
				"Three rows are written to outFile: mean, lower limit, and upper limit",
				"for each of the selected columns.",
				"",
				"Parameters:",
				"	int firstColumn = 0;",
				"	int lastColumn = numeric_limits<int>::max();",
				"	int bootstrap = 200;",
				"	double lowerLimit = 0.05;",
				"	double upperLimit = 0.95;",
				"	int skipLeadingRows = 1;",
				"	int skipTrailingRows = 1;",
				"string inFile, outFile;",
				"",
				"By default, we skip the first row (assumed to contain headings)",
				"and the last row (assumed to contain averages or some such)."
			};

			for (auto & h : help) {
				cerr << h << '\n';
			}
		}

		Params(Args & args) {
			bool ok = true;

			if (args.IsDefined("help")) {
				Help();
				return;
			}

			if (!args.Get("firstColumn", firstColumn)) {
				cerr << "Information: Argument '--firstColumn' not supplied. Using default value.\n";
			}

			if (!args.Get("lastColumn", lastColumn)) {
				cerr << "Information: Argument '--lastColumn' not supplied. Using default value.\n";
			}

			if (lastColumn < firstColumn) {
				lastColumn = firstColumn;
			}

			if (!args.Get("bootstrap", bootstrap)) {
				cerr << "Information: Argument '--bootstrap' not supplied. Using default value.\n";
			}

			if (!args.Get("skipLeadingRows", skipLeadingRows)) {
				cerr << "Information: Argument '--skipLeadingRows' not supplied. Using default value.\n";
			}

			if (skipLeadingRows < 0) {
				skipLeadingRows = 0;
			}

			if (!args.Get("skipTrailingRows", skipTrailingRows)) {
				cerr << "Information: Argument '--skipTrailingRows' not supplied. Using default value.\n";
			}

			if (skipTrailingRows < 0) {
				skipTrailingRows = 0;
			}

			if (!args.Get("lowerLimit", lowerLimit)) {
				cerr << "Information: Argument '--lowerLimit' not supplied. Using default value.\n";
			}

			if (lowerLimit < 0) {
				lowerLimit = 0;
			}

			if (!args.Get("upperLimit", upperLimit)) {
				cerr << "Information: Argument '--upperLimit' not supplied. Using default value.\n";
			}

			if (upperLimit < 0) {
				upperLimit = 0;
			}

			if (!args.Get("inFile", inFile)) {
				cerr << "Argument '--inFile' not supplied.\n";
				ok = false;
			}

			if (!args.Get("outFile", outFile)) {
				cerr << "Argument '--outFile' not supplied.\n";
				ok = false;
			}

			if (!ok) {
				Help();
				throw Exception("Invalid arguments.", FileAndLine);
			}
		}
	};

	struct ColumnAverage {
		static void Run(Args & args) {
			Params parms(args);
			vector<vector<string>> rows;

			ifstream inFile(parms.inFile);
			CsvReader reader(inFile, '\t');
			reader.Read(rows);

			uint lastColumn = std::min(rows[0].size() - 1, (size_t)parms.lastColumn);
			uint wantedCols = lastColumn - parms.firstColumn + 1;

			if (wantedCols <= 0) {
				throw Exception("The number of columns selected is less than 1.\n", FileAndLine);
			}

			uint wantedRows = rows.size() - parms.skipLeadingRows - parms.skipTrailingRows;

			if (wantedRows <= 0) {
				throw Exception("The number of rows selected is less than 1.\n", FileAndLine);
			}

			vector<vector<double>> cols(wantedCols);

			for (uint col = parms.firstColumn; col <= lastColumn; col++) {
				cols[col - parms.firstColumn].resize(wantedRows);

				for (uint i = 0; i < wantedRows; i++) {
					istringstream str(rows[i + parms.skipLeadingRows][col]);
					double value;
					str >> value;
					cols[col - parms.firstColumn][i] = value;
				}
			}

			vector<double> overallMean(cols.size());

			for (uint col = 0; col < cols.size(); col++) {
				overallMean[col] = 0;

				for (uint i = 0; i < wantedRows; i++) {
					overallMean[col] += cols[col][i];
				}

				overallMean[col] /= wantedRows;
			}

			vector<vector<double>> sampleMeans(cols.size());
			UniformIntRandom<size_t> rand(time(0), 0, wantedRows-1);

			for (uint col = 0; col < cols.size(); col++) {
				for (uint repeat = 0; repeat < parms.bootstrap; repeat++) {
					double sum = 0;

					for (uint i = 0; i < wantedRows; i++) {
						size_t sampleIdx = rand();
						sum += cols[col][sampleIdx];
					}

					sampleMeans[col].push_back(sum / wantedRows);
				}
			}

			vector<double> lowerLimits(cols.size());
			vector<double> upperLimits(cols.size());

			for (uint col = 0; col < cols.size(); col++) {
				std::sort(sampleMeans[col].begin(), sampleMeans[col].end());
				lowerLimits[col] = sampleMeans[col][parms.lowerLimit * parms.bootstrap];
				upperLimits[col] = sampleMeans[col][parms.upperLimit * parms.bootstrap];
			}

			ofstream outFile(parms.outFile);
			outFile << "Mean\tLowerLimit\tUpper Limit\n";

			for (uint i = 0; i < overallMean.size(); i++) {
				outFile << overallMean[i] << "\t" << lowerLimits[i] << "\t" << upperLimits[i] << "\n";
			}
		}
	};


};

int main(int argc, char** argv) {
	try {
		Args args(argc, argv);
		Bootstrap::ColumnAverage::Run(args);
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


