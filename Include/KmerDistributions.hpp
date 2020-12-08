#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

namespace QutBio {
	#include "Histogram.hpp"
	#include "SimilarityMatrix.hpp"

	class KmerDistributions {
	public:
		static Histogram<double> GetOneMerDistanceDistribution(
			const SimilarityMatrix & blosum62,
			const Histogram<Symbol> & symbolDist,
			Histogram<double> & oneMerDistances
			) {
			oneMerDistances.GetOneMerHistogram<Symbol, Symbol>(
				symbolDist,
				[blosum62]( Symbol x, Symbol y) { return blosum62.MaxValue() - blosum62.Similarity(x, y); }
			);

			return oneMerDistances;
		}

		static Histogram<double> GetOneMerSimilarityDistribution(
			const SimilarityMatrix & blosum62,
			const Histogram<Symbol> & symbolDist,
			Histogram<double> & oneMerDistances
			) {
			oneMerDistances.GetOneMerHistogram<Symbol, Symbol>(
				symbolDist,
				[blosum62]( Symbol x, Symbol y) { return blosum62.Similarity(x, y); }
			);

			return oneMerDistances;
		}

		static void GetHousdorffAverageFragmentDistributions(
			int maxK,
			int fragLength,
			const Histogram<double> & oneMerDistances,
			map<int, DiscreteDistribution> & hausdorffFragmentDistributions
			) {
			Histogram<double> kmerDistances = oneMerDistances;

			for ( int k = 2; k <= maxK; k++ ) {
				Histogram<double> newHistogram;
				kmerDistances.DoConvolution(oneMerDistances, newHistogram);
				kmerDistances = newHistogram;
				DiscreteDistribution discrete;
				DiscreteDistribution minDist;
				discrete.SetPmf(kmerDistances);
				discrete.GetMinimumDistribution(fragLength, minDist);

				Histogram<double> currentSum = minDist.Pmf();

				for ( int f = 2; f <= fragLength; f++ ) {
					Histogram<double> newSum;
					currentSum.DoConvolution(minDist.Pmf(), newSum);
					newSum.Cleanup([](double key, double value) { return value <= 0; });
					currentSum = newSum;
				}

				Histogram<double> averagePmf;

				for ( auto p : currentSum.data ) {
					averagePmf.data[p.first / fragLength] = p.second;
				}

				DiscreteDistribution averageDistribution;
				averageDistribution.SetPmf(averagePmf);

				DiscreteDistribution hausdorffAverageDistribution;
				averageDistribution.GetMaximumDistribution(2, hausdorffAverageDistribution);
				hausdorffAverageDistribution.Cleanup();
				hausdorffFragmentDistributions[k] = hausdorffAverageDistribution;
			}
		}

		static void GetHausdorffAverageFragmentDistributions(
			int maxK,
			int fragLength,
			const SimilarityMatrix & similarityMatrix,
			const Histogram<Symbol> & symbolDist,
			map<int, DiscreteDistribution> & hausdorffFragmentDistributions
			) {
			Histogram<double> oneMerDistances;
			GetOneMerDistanceDistribution(similarityMatrix, symbolDist, oneMerDistances);
			GetHousdorffAverageFragmentDistributions(maxK, fragLength, oneMerDistances, hausdorffFragmentDistributions);
		}
	};
}