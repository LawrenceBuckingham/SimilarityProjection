#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cstring>

#include "Exception.hpp"

namespace QutBio {
	template<typename T>
	class Subrange {
		T min, max;

		struct ForwardIterator {
			T current;

			ForwardIterator( T current ) : current( current ) {}

			T operator * () { return current; }

			void operator ++ () { current++; }

			void operator ++ (int) { current++; }

			friend bool operator == ( const ForwardIterator & lhs,  const ForwardIterator & rhs ) {
				return lhs.current == rhs.current;
			}

			friend bool operator != ( const ForwardIterator & lhs,  const ForwardIterator & rhs ) {
				return lhs.current != rhs.current;
			}
		};
	public:
		Subrange( const T& min, const T& max ) : min( min ), max( max ) {
			if ( max < min ) {
				throw ArgumentException( FileAndLine );
			}
		}

		ForwardIterator begin() const {
			return ForwardIterator(min);
		}

		ForwardIterator end() const {
			T res = max;
			res++;
			return ForwardIterator(res);
		}

		T Min() const {
			return min;
		}

		T Max() const {
			return max;
		}
	};
}