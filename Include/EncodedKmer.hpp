#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <climits>
#include <cstdint>
#include "Array.hpp"

namespace QutBio {
	typedef uint64_t KmerWord;
	typedef KmerWord * EncodedKmer;
	#define KWORD_BITS (sizeof(KmerWord)*CHAR_BIT)
}
