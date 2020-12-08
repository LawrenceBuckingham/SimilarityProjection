#pragma once

#include "EncodedKmer.hpp"
#include "DistanceType.hpp"

namespace QutBio {

	// Making this an object instead of a static function will let me
	// use it in the KmerClusterAL template class.
	// 2019-05-26: Convert to sparse-coded rather than packed encoding.
	class DnaDistance {
#if INSTRUMENT_DNA_DIST
		size_t callCounter = 0;
#endif
	public:

		// DNA kmers are always encoded in a single 64 bit word,
		// so we can just calculate the Hamming distance by bit 
		// operations.
		Distance GetDistance1(KmerWord x, KmerWord y) const {
			return 64 - POPCOUNT(x & y);
		}

		Distance operator()(KmerWord *x, KmerWord *y, uint kmerLength) const {
#if INSTRUMENT_DNA_DIST
			callCounter++;
#endif
			Distance dist = kmerLength;

			for (size_t i = 0; i * (CHAR_BIT * sizeof(KmerWord) / 4) < kmerLength; i++) {
				uint64_t a = x[i] & y[i];
				dist -= POPCOUNT(a);
			}

			return dist;
		}

		/*
		**	Compute the distance between (presumed lowercase) kmers
		**	via direct comparison.
		*/
		Distance operator() (const char *x, const char * y, uint kmerLength) {
			uint sim = 0;

			for (uint i = 0; i < kmerLength; i++ ) {
				char xc = x[i];
				char yc = y[i];

				if (xc != 'n' && yc != 'n') {
					sim += xc == yc;
				}
			}

			return kmerLength - sim;
		}

		/*
		**	Compute the distance between kmers encoded using sparse code.
		**	sim(\alpha,\beta) = \alpha & beta;
		*/
		Distance operator() (const uint8_t *x, const uint8_t * y, uint kmerLength) {
			uint sim = 0;

			for (uint i = 0; i < kmerLength; i++ ) {
				char xc = x[i];
				char yc = y[i];

				sim += ((xc & yc) != 0);
			}

			return kmerLength - sim;
		}

#if INSTRUMENT_DNA_DIST
		size_t CallCounter() {
			return callCounter;
		}
#endif

	};
}
