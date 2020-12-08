#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <atomic>
#include <omp.h>

#include "Array.hpp"
#include "FastaSequence.hpp"
#include "Fragment.hpp"
#include "SimilarityMatrix.hpp"

namespace QutBio {

	template<typename vectorType, typename distanceType>
	class FragDistProcessor {
	private:
		omp_lock_t lock;

	protected:
		vectorType & empty() {
			static vectorType instance;
			return instance;
		}

	public:
		vectorType rowMinima;
		vectorType colMinima;
		size_t queryKmerCount;
		size_t subjectKmerCount;

		FragDistProcessor(
			vectorType &rowMinima,
			vectorType &colMinima,
			size_t queryKmerCount,
			size_t subjectKmerCount
		) : rowMinima(rowMinima),
			colMinima(colMinima),
			queryKmerCount(queryKmerCount),
			subjectKmerCount(subjectKmerCount)
			//
		{
			omp_init_lock(&lock);
		}

		virtual ~FragDistProcessor() {
			omp_destroy_lock(&lock);
		}

		virtual void ProcessDistance(size_t queryPos, size_t subjectPos, Distance distance) = 0;

		void Reset(size_t queryKmerCount, Distance defaultDistance) {
			this->queryKmerCount = queryKmerCount;
			size_t n = colMinima.size();
			std::fill_n(colMinima.data(), n, defaultDistance);
			std::fill_n(rowMinima.data(), queryKmerCount, defaultDistance);
		}

		omp_lock_t * Lock() {
			return &lock;
		}
	};

	template<typename vectorType, typename distanceType>
	class FragDistProcessor_Frag : public FragDistProcessor<vectorType, distanceType> {
	public:
		size_t queryFragCount;
		size_t subjectFragCount;
		size_t fragmentLength;
		size_t fragsPerTile;

		int hits;

		FragDistProcessor_Frag() : FragDistProcessor<vectorType, distanceType>(this->empty(), this->empty(), 0, 0) {}

		FragDistProcessor_Frag(
			vectorType &rowMinima,
			vectorType &colMinima,
			size_t queryKmerCount,
			size_t subjectKmerCount,
			size_t queryFragCount,
			size_t subjectFragCount,
			size_t fragmentLength,
			int	fragsPerTile,
			size_t *queryFragMapping,
			size_t *subjectFragMapping
		) :
			FragDistProcessor<vectorType, distanceType>(rowMinima, colMinima, queryKmerCount, subjectKmerCount),
			queryFragCount(queryFragCount),
			subjectFragCount(subjectFragCount),
			fragmentLength(fragmentLength),
			fragsPerTile(fragsPerTile) {}

		void ProcessDistance(size_t queryPos, size_t subjectPos, Distance distance) {
			omp_set_lock(FragDistProcessor<vectorType, distanceType>::Lock());

			hits++;

			size_t iMax = queryPos / fragmentLength;
			size_t jMax = subjectPos / fragmentLength;

			//if ( iMax >= data->queryFragCount || jMax >= data->subjectFragCount ) {
			//	cerr << "kmer is not included in a fragment!" << endl;
			//	throw Exception("This is not supposed to happen.", FileAndLine);
			//}

			if (fragsPerTile == 1) {
				if (distance < this->rowMinima[iMax]) this->rowMinima[iMax] = distance;
				if (distance < this->colMinima[jMax]) this->colMinima[jMax] = distance;
			}
			else {
				cerr << "Non-unit fragsPerTile is temporarily out of service. Sorry.\n";
				exit(1);

				for ( size_t i = 0; i <= iMax && i < fragsPerTile; i++) {
					if (distance < this->rowMinima[iMax - i]) {
						this->rowMinima[iMax - i] = distance;
					}
				}
				for ( size_t j = 0; j <= jMax && j < fragsPerTile; j++) {
					if (distance < this->colMinima[jMax - j]) {
						this->colMinima[jMax - j] = distance;
					}
				}
			}

			omp_unset_lock(FragDistProcessor<vectorType, distanceType>::Lock());
		}

		/**
		 *	<summary>
		 *		Populates the array referred to by fragMapping, which must have kmerCount elements allocated,
		 *		with the index of the right-most partition to which each of the kmers in the sequence belongs.
		 *		If fragsPerTile is greater than 1, then kmer i belongs to fragments fragMapping[i], fragMapping[i-1],
		 *		..., fragMapping[i-fragsPerTile+1].
		 *	</summary>
		 *	<param name="fragMapping">An array that will be populated with the fragment index values, as noted above.</param>
		 *	<param name="kmerCount">The number of kmers in the sequence.</param>
		 *	<param name="fragCount">The number of fragments required to partition the sequence.</param>
		 *	<param name="stepSize">The equalised fragment length.</param>
		 */

		static void GetFragmentMapping(vector<size_t> & fragMapping, size_t kmerCount, size_t fragCount, double stepSize) {
			uint start = 0;

			for (int f = 0; f < fragCount; f++) {
				uint end = Fragment::GetFragmentStart(f + 1, stepSize, kmerCount);

				for (uint j = start; j < end; j++) {
					fragMapping[j] = f;
				}

				start = end;
			}
		}
	};

	template<typename vectorType, typename distanceType>
	class FragDistProcessor_NoFrag : public FragDistProcessor<vectorType, distanceType> {
	public:
		size_t skip = 0;

		FragDistProcessor_NoFrag() : FragDistProcessor<vectorType, distanceType>(this->empty(), this->empty(), 0, 0) {}

		FragDistProcessor_NoFrag(
			vectorType &rowMinima,
			vectorType &colMinima,
			size_t queryKmerCount,
			size_t subjectKmerCount,
			size_t skip
		) :
			FragDistProcessor<vectorType, distanceType>(rowMinima, colMinima, queryKmerCount, subjectKmerCount), skip(skip) {
			Assert::IsTrue(skip > 0, FileAndLine);
		}

		void ProcessDistance(size_t queryPos, size_t subjectPos, Distance distance) {
			omp_set_lock(FragDistProcessor<vectorType, distanceType>::Lock());
			size_t queryPos2 = queryPos / skip;
			if (distance < this->rowMinima[queryPos2]) this->rowMinima[queryPos2] = distance;
			if (distance < this->colMinima[subjectPos]) this->colMinima[subjectPos] = distance;
			omp_unset_lock(FragDistProcessor<vectorType, distanceType>::Lock());
		}
	};

}
