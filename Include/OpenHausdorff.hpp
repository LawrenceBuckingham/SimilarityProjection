#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "FragmentAggregationMode.hpp"
#include "HausdorffCalculator.hpp"

namespace QutBio {
	/**
	 *	<summary>
	 *		A subclass of HausdorffCalculator that reveals the mutual distance matrix and fragment distance cache
	 *		for use by programs that need efficient access to an exhaustive MxN distance matrix.
	 */
	class OpenHausdorff :
		public HausdorffCalculator {
	public:
		OpenHausdorff(
			SimilarityMatrix * matrix,
			uint kmerLength,
			FragmentAggregationMode * kmerMode,
			FragmentAggregationMode * fragMode,
			Alphabet * alphabet,
			uint fragLength,
			uint maxQueryLength,
			uint maxSubjectLength
			) : HausdorffCalculator(
			matrix, kmerLength, kmerMode, fragMode, alphabet, fragLength, maxQueryLength, maxSubjectLength
			) {}

		const MatrixView<Distance> & Distances() const { return kmerDistCache; }

		vector<Distance> & RowMinima() { return rowMinima; }
		vector<Distance> & ColMinima() { return colMinima; }
		FragmentDistance FragmentDistance() { return fragmentDistance; }
		CollectionDistance CollectionDistance () { return collectionDistance; }
	};

}
