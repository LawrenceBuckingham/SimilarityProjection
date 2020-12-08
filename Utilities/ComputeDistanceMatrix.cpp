#include <map>
#include <string>
#include <vector>
#include <algorithm>

#include <Args.hpp>
#include <FastaSequence.hpp>
#include <OmpTimer.h>
#include <Exception.hpp>
#include <FileUtil.hpp>
#include <HausdorffCalculator.hpp>
#include <DataLoader.hpp>

// file://a:\phd\Notes\DynamicProgrammingAlgorithm.html

using namespace std;
using namespace QutBio;

namespace Utilities {
	struct ComputeDistanceMatrix {

		struct ReduceMode : public EnumBase {
		private:
			ReduceMode(string literal, int value) : EnumBase(literal, value) {}
			static std::mutex & m() {
				static std::mutex mut;
				return mut;
			}

		public:
			static ReduceMode * Hausdorff() {
				std::unique_lock < mutex > lck{ m() };
				static ReduceMode value("hausdorff", 0);
				return &value;
			}

			static ReduceMode * Jaccard() {
				std::unique_lock < mutex > lck{ m() };
				static ReduceMode value("jaccard", 1);
				return &value;
			}

			static vector<EnumBase *> Values() {
				std::vector<EnumBase *> result{
					Hausdorff(),
					Jaccard()
				};
				return result;
			}
		};

		struct Params {
			bool ok;
			string query;
			string db;
			string dist;
			string kmerDist;
			size_t idIndex;
			size_t kmerLength;
			size_t maxRecords = 500;
			Distance threshold;
			Alphabet * alphabet;
			SimilarityMatrix * similarity;
			ReduceMode *mode;

			static void Help() {
				vector<string> help = {
"ComputeDistanceMatrix: calculates the pairwise similarity, plus optional k-mer distance matrices ",
"Arguments:",
"--help: display this text.",
"--matrix file-name:"
"--threshold t:"
"\n	a numeric value (integer valued in the present version) that is "
"\n	used to classify k-mer pairs as being sufficiently close to count "
"\n	as instances of the same word: d(x,y) <= t ==> treat-as-same.",
"--mode hausdorff | jaccard:"
"\n	the reduction mode used to summarise the distance matrix.",
"--maxRecords int:"
"\n	Optional. The maximum number of ranked hits to output. Default"
"\n	value is 500."
				};
				for (auto & s : help) {
					cerr << s << "\n\n";
				}
			}

			Params(Args & args) {

#define REQ(x,h) args.Required(x, #x, h)
#define OPT(x,h) args.Optional(x, #x, h)
#define STR(x) STR2(x)
#define STR2(x) #x

				REQ(query, "The file path of a FASTA formatted query dataset.");

				REQ(db, "The file path of a FASTA formatted reference dataset.");

				REQ(dist,
					"The file path of a document that will be populated with the pairwise distance \n"
					"matrix. Each entry on the matrix is a tuple containing (query_id, subject_id, \n"
					"distance).");

				OPT(kmerDist,
					"The file path of a document that will be populated with a pairwise k-mer \n"
					"distance matrix for each (query,subject) pair. Records are written over \n"
					"multiple lines. The first line of a record contains a record of the form: \n"
					"query-id  subject-id  distance. This is followed by the matrix, written as \n"
					"(m+1-k) rows, each containing (n+1-k) numeric values separated by tab \n"
					"characters. Here, m is the character length of the query sequence and n is \n"
					"the character length of the reference sequence. The matrix is followed by a \n"
					"single empty line. This generates a huge volume of data, so use it only on \n"
					"small databases.");

				OPT(maxRecords, "The maximum number of ranked matches to return.");

				REQ(kmerLength, "The word length for k-mer tiling.");

				REQ(idIndex,
					"The zero-origin location of the sequence Id within the pipe-separated list of \n"
					"metadata fields in the definition line.");

				args.Required("mode", ReduceMode::Values(), mode,
					"The distance calculation for kmer comparison. (hausdorff|jaccard)");

				if (mode == ReduceMode::Jaccard()) {
					REQ(threshold,
						"The threshold distance for determining if two k-mers are to be considered to \n"
						"be instances of a common vocabulary term.");
				}

				args.Required(alphabet,similarity);

				if (!args.Ok()) {
					args.Help();
					throw Exception("Invalid arguments.", FileAndLine);
				}
			}
		};

		struct Ranking {
			const char * id;
			double distance;

			Ranking(
				const char * id = nullptr, double distance = 0.0
			) :
				id(id), distance(distance) {}
		};

		static void Run(int argc, char** argv) {
			Args args(argc, argv);
			Params parms(args);
			Alphabet &alphabet = *parms.alphabet;

			auto db = Load::Fasta(parms.db, parms.idIndex, &alphabet);
			PadSequences(db, parms.kmerLength, alphabet.DefaultSymbol());

			auto query = Load::Fasta(parms.query, parms.idIndex, &alphabet);
			PadSequences(db, parms.kmerLength, alphabet.DefaultSymbol());

			size_t maxQueryLen = 0;

			for (auto seq : query) {
				size_t n = seq->KmerCount(parms.kmerLength);

				if (n > maxQueryLen) {
					maxQueryLen = n;
				}
			}

			size_t maxDbLen = 0;

			for (auto seq : db) {
				size_t n = seq->KmerCount(parms.kmerLength);

				if (n > maxDbLen) {
					maxDbLen = n;
				}
			}

			Distance distLookup[128][128];
			parms.similarity->PopulateDistanceTable(distLookup);

			ofstream distFile(parms.dist);
			ofstream *matrixFile = parms.kmerDist.length() > 0 ? new ofstream(parms.kmerDist) : nullptr;

#pragma omp parallel
			{
				vector<Distance> distMatrix(maxQueryLen * maxDbLen);
				MatrixView<Distance> matrix(distMatrix.data());
				vector<Ranking> distances(db.size());
				vector<Distance> rowSummary(maxQueryLen);
				vector<Distance> colSummary(maxDbLen);

#pragma omp for
				for (size_t queryPos = 0; queryPos < query.size(); queryPos++) {
					auto & querySeq = *query[queryPos];
					size_t m = querySeq.KmerCount(parms.kmerLength);
					distances.clear();

					for (size_t dbPos = 0; dbPos < db.size(); dbPos++) {
						auto & dbSeq = *db[dbPos];
						size_t n = dbSeq.KmerCount(parms.kmerLength);
						matrix.Reinterpret(m, n);
						HausdorffCalculator::ComputeDistanceMatrix(querySeq.Sequence().data(), dbSeq.Sequence().data(), m, n, parms.kmerLength, distLookup, matrix);
						double distance = parms.mode == ReduceMode::Jaccard()
							? ComputeDistanceJaccard(matrix, rowSummary, colSummary, parms.threshold)
							: ComputeDistanceHausdorff(matrix, rowSummary, colSummary);

#if 0
						if (querySeq.IdStr() == "g00000" && dbSeq.IdStr() == "g17769") {
#pragma omp critical
							{
								cerr << querySeq.IdStr() << " " << dbSeq.IdStr() << "\n";
								cerr << "rowMinima";
								for (uint i = 0; i < m; i++) { auto d = rowSummary[i]; cerr << "\t" << d; }
								cerr << "\n";
								cerr << "colMinima";
								for (uint i = 0; i < n; i++) { auto d = colSummary[i]; cerr << "\t" << d; }
								cerr << "\n";

								abort();
							}
						}
#endif

						distances.emplace_back(dbSeq.IdStr().c_str(), distance);

						if (matrixFile) {
#pragma omp critical 
							{
								(*matrixFile) << querySeq.IdStr() << "\t" << dbSeq.IdStr() << "\t" << distance << "\n";
								for (size_t r = 0; r < m; r++) {
									(*matrixFile) << matrix(r, 0);
									for (size_t c = 1; c < n; c++) {
										(*matrixFile) << "\t" << matrix(r, c);
									}
									(*matrixFile) << "\n";
								}
								(*matrixFile) << "\n";
							}
						}
					}

#pragma omp critical
					{
						std::sort(distances.begin(), distances.end(), [](const auto & lhs, const auto & rhs) { return lhs.distance < rhs.distance; });
						size_t i = 0;
						for (auto & d : distances) {
							if (i >= parms.maxRecords) break;

							distFile << querySeq.IdStr() << "\t" << d.id << "\t" << d.distance << "\n";
							i++;
						}
					}
				}
			}

			if (matrixFile) delete matrixFile;
		}

		static void PadSequences(
			vector<FastaSequence *> & seqs,
			size_t minLength,
			Symbol padding
		) {
			for (auto seq : seqs) {
				seq->Pad(minLength, padding);
			}
		}

		static double ComputeDistanceHausdorff(
			MatrixView<Distance> &matrix,
			vector<Distance> & rowSummary,
			vector<Distance> & colSummary
		) {
			size_t
				m = matrix.Rows(),
				n = matrix.Cols();

			if (m == 0 || n == 0) return numeric_limits<double>::max();

			for (size_t r = 0; r < m; r++) rowSummary[r] = numeric_limits<Distance>::max();
			for (size_t c = 0; c < n; c++) colSummary[c] = numeric_limits<Distance>::max();

			for (size_t r = 0; r < m; r++) {
				for (size_t c = 0; c < n; c++) {
					auto d = matrix(r, c);
					if (d < rowSummary[r]) rowSummary[r] = d;
					if (d < colSummary[c]) colSummary[c] = d;
				}
			}
			double rowTot = 0, colTot = 0;
			for (size_t r = 0; r < m; r++) rowTot += rowSummary[r];
			for (size_t c = 0; c < n; c++) colTot += colSummary[c];
			return (rowTot / m + colTot / n) / 2;
		}

		static double ComputeDistanceJaccard(MatrixView<Distance> &matrix, vector<Distance> & rowSummary, vector<Distance> & colSummary, Distance threshold) {
			size_t
				m = matrix.Rows(),
				n = matrix.Cols();

			if (m == 0 || n == 0) return numeric_limits<double>::max();

			for (size_t r = 0; r < m; r++) {
				for (size_t c = 0; c < n; c++) {
					auto d = matrix(r, c);
					if (d < threshold) {
						rowSummary[r] = 1;
						colSummary[c] = 1;
					}
				}
			}

			double rowTot = 0, colTot = 0;
			for (size_t r = 0; r < m; r++) rowTot += rowSummary[r];
			for (size_t c = 0; c < n; c++) colTot += colSummary[c];
			return 1 - (rowTot / m + colTot / n) / 2;
		}
	};
}

using namespace Utilities;

std::mutex DistanceType::m;

int main(int argc, char** argv) {
	try {
		ComputeDistanceMatrix::Run(argc, argv);
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

