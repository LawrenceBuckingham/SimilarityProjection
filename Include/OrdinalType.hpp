#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

namespace QutBio {

	/// <summary>
	/// Implement a type-safe "ordinal" type which is harder to mix up with arithmetic types!
	/// </summary>
	/// <typeparam name="T">The encapsulated numeric type.</typeparam>
	template<typename T>
	class Ordinal {
	public:
		T value;

		Ordinal( T val = T{} ) : value( val ) {}

		static Ordinal From( T value ) {
			T res{ value };
			return res;
		}

		friend bool operator == ( const Ordinal& lhs, const Ordinal& rhs ) {
			return lhs.value == rhs.value;
		}

		friend bool operator != ( const Ordinal& lhs, const Ordinal& rhs ) {
			return lhs.value != rhs.value;
		}

		friend bool operator < ( const Ordinal& lhs, const Ordinal& rhs ) {
			return lhs.value < rhs.value;
		}

		friend bool operator <= ( const Ordinal& lhs, const Ordinal& rhs ) {
			return lhs.value <= rhs.value;
		}

		friend bool operator > ( const Ordinal& lhs, const Ordinal& rhs ) {
			return lhs.value > rhs.value;
		}

		friend bool operator >= ( const Ordinal& lhs, const Ordinal& rhs ) {
			return lhs.value >= rhs.value;
		}

		void operator++() {
			value++;
		}

		void operator++(int) {
			value++;
		}

		template<typename U>
		friend Ordinal operator+( const Ordinal& lhs, const U & delta ) {
			Ordinal res{ (T) ( lhs.value + delta ) };
			return res;
		}

		Ordinal& operator--() {
			value--;
			return *this;
		}

		template<typename U>
		friend Ordinal operator-( const Ordinal& lhs, const U & delta ) {
			Ordinal res{ (T) ( lhs.value - delta ) };
			return res;
		}

		friend std::ostream& operator << ( std::ostream& str, const Ordinal& ch ) {
			return str << ch.value;
		}
	};
}
