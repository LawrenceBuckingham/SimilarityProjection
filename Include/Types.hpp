#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "OrdinalType.hpp"

#include <cstring>

namespace QutBio {
	enum YesNoCancel {
		Cancel, Yes, No
	};

	using byte = unsigned char;
	using uint = unsigned int;
	using ulong = unsigned long;

	using Symbol = Ordinal<byte>;

	// TODO: convert this to Ordinal<unsigned short>
	using Digram = unsigned short;

	using Trigram = unsigned int;


	template<typename T, size_t Capacity = 256>
	class ByteIdxArray {
		T values[Capacity];
	public:
		ByteIdxArray() {
			Clear();
		}

		void Clear() {
			std::memset( values, 0, sizeof( values ) );
		}

		T& operator[]( const Symbol & x ) {
			return values[x.value];
		}

		T at( const Symbol & x ) const {
			return values[x.value];
		}

		T* begin() const { return values; }
		T* end() const { return values + Capacity; }
	};

	template <typename T1, typename T2, typename T3>
	struct Triple {
		T1 item1;
		T2 item2;
		T3 item3;

		Triple() {}

		Triple( const T1& item1, const T2& item2, const T3& item3 ) :
			item1( item1 ), item2( item2 ), item3( item3 ) {}

		friend bool operator<( const Triple& lhs, const Triple& rhs ) {
			if ( lhs.item1 < rhs.item1 ) return true;
			if ( rhs.item1 < lhs.item1 ) return false;
			if ( lhs.item2 < rhs.item2 ) return true;
			if ( rhs.item2 < lhs.item2 ) return false;
			if ( lhs.item3 < rhs.item3 ) return true;
			if ( rhs.item3 < lhs.item3 ) return false;
			return false;
		}

		friend bool operator==( const Triple& lhs, const Triple& rhs ) {
			return lhs.item1 == rhs.item1 && lhs.item2 == rhs.item2 && lhs.item3 == rhs.item3;
		}

		friend bool operator!=( const Triple& lhs, const Triple& rhs ) {
			return !operator==( lhs, rhs );
		}
	};

	template <typename T1, typename T2>
	struct Pair {
		T1 item1;
		T2 item2;

		Pair() {}

		Pair( const T1& item1, const T2& item2 ) :
			item1( item1 ), item2( item2 ) {}

		friend bool operator<( const Pair& lhs, const Pair& rhs ) {
			if ( lhs.item1 < rhs.item1 ) return true;
			if ( rhs.item1 < lhs.item1 ) return false;
			if ( lhs.item2 < rhs.item2 ) return true;
			if ( rhs.item2 < lhs.item2 ) return false;
			return false;
		}

		friend bool operator==( const Pair& lhs, const Pair& rhs ) {
			return lhs.item1 == rhs.item1 && lhs.item2 == rhs.item2;
		}

		friend bool operator!=( const Pair& lhs, const Pair& rhs ) {
			return !operator==( lhs, rhs );
		}
	};
}
