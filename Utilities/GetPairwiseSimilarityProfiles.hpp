#pragma once

#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>

#include "QutBio.h"
#include "Homologs.hpp"
#include "OpenHausdorff.hpp"

using namespace QutBio;

namespace Utilities {

	class GetPairwiseSimilarityProfiles {
	private:
		struct Parameters {
			/// <summary>Name of file containing FASTA-formatted sequences.</summary>
			string dbFile;

			/// <summary>Name of file containing ids of sequences to be processed in this run.</summary>
			string idFile;

			/// <summary>Optional name of file containing custom similarity matrix.</summary>
			string matrixFile;

			/// <summary>0-origin index of sequence id in pipe-separated fasta defline.</summary>
			int idIndex = 0;

			/// <summary>0-origin index of field containing class labels (semicolon-separated if more than 1) 
			///	in pipe-separated fasta defline.
			/// </summary>
			int classIndex = -1;

			/// <summary>K in k-mer.</summary>
			int kmerLength = 0;

			/// <summary>Enumerated value if using one of the built-in similarity matrix types.</summary>
			DistanceType * dist = DistanceType::HalperinEtAl();

			/// <summary>Numeric matrix type if using built-in type.</summary>
			int matrixId = 62;

			/// <summary>True iff the alphabet is case-sensitive.</summary>
			bool isCaseSensitive = false;

			/// <summary>Writes a text representation of the parameters to a stream.</summary>
			friend ostream & operator<<(ostream& out, const Parameters & p) {
				out << "Parameters:" << endl;
				out << "--dbFile " << p.dbFile << endl;
				out << "--idFile " << p.idFile << endl;
				out << "--matrixFile " << p.matrixFile << endl;
				out << "--idIndex " << p.idIndex << endl;
				out << "--classIndex " << p.classIndex << endl;
				out << "--kmerLength " << p.kmerLength << endl;
				out << "--dist " << p.dist->Name() << endl;
				out << "--matrixId " << p.matrixId << endl;
				out << "--isCaseSensitive " << p.isCaseSensitive << endl;
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

			auto alphabet = Alphabets::AA();
			SimilarityMatrix * matrix = SimilarityMatrix::GetMatrix(alphabet, p.dist, p.matrixId, p.matrixFile);

			auto db = Load::Fasta(p.dbFile, p.idIndex, alphabet);
			cerr << args.ProgName() << ": " << db.size() << " sequences loaded.\n";

			Index<FastaSequence> idx(db);
			vector<FastaSequence *> items;

			File::ReadStrings(p.idFile, [&](const string & s) {
				auto seqIter = idx.find(s);

				if (seqIter != idx.end()) {
					items.push_back(seqIter->second);
				}
			});

			const size_t K = p.kmerLength;

			size_t maxKmerCount = accumulate(items.begin(), items.end(), size_t(0), [=](size_t maxLen, FastaSequence * seq) {
				return std::max(maxLen, seq->KmerCount(K));
			});

#pragma omp parallel for
			for (int i = 0; i < (int)items.size(); i++) {
				auto & query = *items[i];
				auto queryChars = query.Sequence().data();
				const size_t M = query.KmerCount(K);

				for (int j = i; j < (int)items.size(); j++) {
					Projector distanceCalc(matrix, p.kmerLength, FragmentAggregationMode::HausdorffAverageAverage());
					vector<Distance> rowMinima(maxKmerCount);
					vector<Distance> colMinima(maxKmerCount);

					auto & subject = *items[j];
					auto subjectChars = subject.Sequence().data();
					const size_t N = subject.KmerCount(K);

					fprintf(stderr, "Processing (%s,%s)\n", query.IdStr().c_str(), subject.IdStr().c_str());

					distanceCalc.ComputeDistanceMatrix( queryChars, subjectChars, rowMinima, colMinima, M, N );

#pragma omp critical
					{
						WriteSimilarity(query.IdStr(), subject.IdStr(), M, N, rowMinima, colMinima);
					}
				}
			}

			vector<const EncodedFastaSequence *> equivalenceClasses;
		}

		static void WriteSimilarity(
			const string & queryId,
			const string & subjectId,
			const size_t M,
			const size_t N,
			const vector<Distance> &rowMinima,
			const vector<Distance> &colMinima
		) {
			cout << "query," << queryId << "\n";
			cout << "subject," << subjectId << "\n";
			cout << "rowMinima"; for (uint i = 0; i < M; i++) { cout << "," << rowMinima[i]; } cout << "\n";
			cout << "colMinima"; for (uint i = 0; i < N; i++) { cout << "," << colMinima[i]; } cout << "\n";
			cout << "\n";
		}

		/// <summary> Parses and validates the arguments, returning the results in a Parameters object.
		/// </summary>
		/// <param name="arguments"></param>
		/// <returns></returns>

		static void GetValidatedParameters(Args & arguments, Parameters & p) {
			if (!arguments.Get("dbFile", p.dbFile)) {
				cerr << "Argument 'dbFile' not defined." << endl;
				abort();
			}

			if (!arguments.Get("idFile", p.idFile)) {
				cerr << "Argument 'idFile' not defined." << endl;
				abort();
			}

			if (arguments.IsDefined("idIndex")) {
				arguments.Get("idIndex", p.idIndex);
			}

			if (!(p.idIndex >= 0)) {
				cerr << "Argument 'idIndex' not valid." << endl;
				abort();
			}

			if (arguments.IsDefined("classIndex")) {
				arguments.Get("classIndex", p.classIndex);
			};

			if (!(p.classIndex != p.idIndex)) {
				cerr << "Argument 'classIndex' must be different from 'idIndex'." << endl;
				abort();
			}

			if (!(arguments.Get("kmerLength", p.kmerLength))) {
				cerr << "Argument 'kmerLength' not valid." << endl;
				abort();
			}

			if (!(p.kmerLength > 0)) {
				cerr << "Argument 'kmerLength' not valid." << endl;
				abort();
			}

			if (!(arguments.Get("dist", DistanceType::Values(), p.dist))) {
				cerr << "Argument 'dist' not valid." << endl;
				abort();
			}

			if (arguments.IsDefined("matrixId")) {
				if (!(arguments.Get("matrixId", p.matrixId))) {
					cerr << "Argument 'matrixId' not valid." << endl;
					abort();
				}
			}
			else {
				p.matrixId = 62;
			}

			if (p.dist == DistanceType::Custom()) {
				if (!arguments.IsDefined("matrixFile")) {
					cerr << "Argument 'matrixFile' not supplied for custom similarity matrix." << endl;
					abort();
				}
				else {
					arguments.Get("matrixFile", p.matrixFile);
				}
			}
			else {
				vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

				bool found = false;

				for (auto x : matrices) {
					if (x == p.matrixId) { found = true; }
				}

				if (!found) {
					cerr << "Matrix id not recognised." << endl;
					abort();
				}
			}

			if (arguments.IsDefined("isCaseSensitive")) {
				if (!(arguments.Get("isCaseSensitive", p.isCaseSensitive))) {
					cerr << "Argument 'isCaseSensitive' not valid." << endl;
					abort();
				}
			}
		}
	};
}
