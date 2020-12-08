#pragma once

#include <numeric>
#include <set>
#include <iomanip>
#include <omp.h>

#include "QutBio.h"
#include "Homologs.hpp"

using namespace QutBio;

namespace Utilities {

	class GetKmerDistanceSample {
	private:
		struct Parameters {
			string matrixFile;
			string outFile;

			int kmerLength = 0;

			DistanceType * dist = DistanceType::HalperinEtAl();
			int matrixId = 62;

			time_t seed = time(0);
			bool isCaseSensitive = false;

			uint sampleSize = 10000;
		};
	public:

		/**
		Processes a ranking file produced by KmerRank to construct intra-class and
		inter-class empirical distributions of distances.
		The format of the input dataset is (laid out in comma-separated columns):

		query_id,ignored_class_id,subject_id,ignored_class_id,distance,rank,isHomolog

		*/
		static void Run(Args args) {
			Parameters parms;
			GetValidatedParameters(args, parms);

			SimilarityMatrix * matrix = SimilarityMatrix::GetMatrix(Alphabets::AA(), parms.dist, parms.matrixId, parms.matrixFile);

			string alphabet = matrix->Alphabet()->Symbols();
			UniformIntRandom<size_t> rand((int) parms.seed, 0, alphabet.size() - 1);

			ofstream outFile(parms.outFile);

			for ( size_t i = 0; i < parms.sampleSize; i++ ) {
				string k1, k2;
				GenerateKmer(parms.kmerLength, alphabet, rand, k1);
				GenerateKmer(parms.kmerLength, alphabet, rand, k2);
				auto dist = matrix->Difference(k1.c_str(), k2.c_str(), parms.kmerLength);
				outFile << k1 << "," << k2 << "," << dist << endl;
			}

			outFile.close();
		}

		static void GenerateKmer(
			size_t kmerLength, 
			const string & alphabet, 
			UniformIntRandom<size_t> &rand, 
			string & kmer
		) {
			kmer.resize(kmerLength);

			for ( uint i = 0; i < kmerLength; i++ ) {
				kmer[i] = alphabet[rand()];
			}
		}

		/// <summary> Parses and validates the arguments, returning the results in a Parameters object.
		/// </summary>
		/// <param name="arguments"></param>
		/// <returns></returns>

		static void GetValidatedParameters(Args & arguments, Parameters & parms) {
			if ( !arguments.Get("outFile", parms.outFile) ) {
				cerr << "Argument 'outFile' not defined." << endl;
				abort();
			}

			if ( !(arguments.Get("kmerLength", parms.kmerLength)) ) {
				cerr << "Argument 'kmerLength' not valid." << endl;
				abort();
			}

			if ( !(parms.kmerLength > 0) ) {
				cerr << "Argument 'kmerLength' not valid." << endl;
				abort();
			}

			if ( !(arguments.Get("sampleSize", parms.sampleSize)) ) {
				cerr << "Argument 'sampleSize' not valid." << endl;
				abort();
			}

			if ( !(parms.sampleSize > 0) ) {
				cerr << "Argument 'sampleSize' must be greater than zero." << endl;
				abort();
			}

			if ( !(arguments.Get("dist", DistanceType::Values(), parms.dist)) ) {
				cerr << "Argument 'dist' not valid." << endl;
				abort();
			}

			if ( arguments.IsDefined("matrixId") ) {
				if ( !(arguments.Get("matrixId", parms.matrixId)) ) {
					cerr << "Argument 'matrixId' not valid." << endl;
					abort();
				}
			}
			else {
				parms.matrixId = 62;
			}

			if ( parms.dist == DistanceType::Custom() ) {
				if ( !arguments.IsDefined("matrixFile") ) {
					cerr << "Argument 'matrixFile' not supplied for custom similarity matrix." << endl;
					abort();
				}
				else {
					arguments.Get("matrixFile", parms.matrixFile);
				}
			}
			else {
				vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

				bool found = false;

				for ( auto x : matrices ) {
					if ( x == parms.matrixId ) { found = true; }
				}

				if ( !found ) {
					cerr << "Matrix id not recognised." << endl;
					abort();
				}
			}

			if ( arguments.IsDefined("seed") ) {
				int seed;

				if ( !(arguments.Get("seed", seed)) ) {
					cerr << "Argument 'seed' not valid." << endl;
					abort();
				}

				if ( seed > 0 ) parms.seed = seed;
			}

			if ( arguments.IsDefined("isCaseSensitive") ) {
				if ( !(arguments.Get("isCaseSensitive", parms.isCaseSensitive)) ) {
					cerr << "Argument 'isCaseSensitive' not valid." << endl;
					abort();
				}
			}
		}
	};
}
