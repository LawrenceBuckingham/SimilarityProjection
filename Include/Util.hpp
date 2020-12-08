#pragma once

//  Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <cstdint>
#include <functional>
#include <cmath>
#include <vector>
#include <unordered_map>

#include "Types.hpp"
#include "Exception.hpp"

#ifndef STRING
#define STRING(x) STRING2(x)
#define STRING2(x) #x
#endif

#define NAMEOF_VARIABLE(Variable) (void(Variable),#Variable)
#define NAMEOF_FUNCTION(Function) (void(&Function),#Function)
#define NAMEOF_METHOD(ClassName,Method) (void(&ClassName::Method),#Method)
#define NAMEOF_TYPE(Type) (void(sizeof(Type)),#Type)

using namespace std;

namespace QutBio {

	//using string = std::string;
	//using istringstream = std::istringstream;
	//using ostringstream = std::ostringstream;
	//using istream = std::istream;
	//using ifstream = std::ifstream;
	//using ostream = std::ostream;

	template <typename T>
	T min( T x, T y ) { return x < y ? x : y; }

	template <typename T>
	T max( T x, T y ) { return x > y ? x : y; }

	template<typename T, typename Iter>
	static T sum( Iter begin, Iter end, T initial ) {
		T total = initial;

		for ( auto iter = begin; iter != end; iter++ ) {
			total += *iter;
		}

		return total;
	}

	template<typename T, typename Collection>
	void GetMin( const Collection& c, T& t ) {
		auto iter = c.begin();
		t = *iter++;

		for ( ; iter != c.end(); iter++ ) {
			if ( *iter < t ) t = *iter;
		}
	}

	template<typename T, typename Collection>
	void GetMax( const Collection& c, T& t ) {
		auto iter = c.begin();
		t = *iter++;

		for ( ; iter != c.end(); iter++ ) {
			if ( *iter > t ) t = *iter;
		}
	}

	class Util {
	public:
		template<typename CollectionType>
		static bool Equal( const CollectionType& lhs, const CollectionType& rhs ) {
			auto lhsBegin = lhs.begin();
			auto lhsEnd = lhs.end();
			auto rhsBegin = rhs.begin();
			auto rhsEnd = rhs.end();
			auto i = lhsBegin;
			auto j = rhsBegin;

			while ( i != lhsEnd && j != rhsEnd ) {
				if ( *i != *j ) return false;
				i++;
				j++;
			}

			return ( i == lhsEnd && j == rhsEnd );
		}

		template<typename T, typename Iter>
		static T Max( Iter begin, Iter end, T initial ) {
			T best = initial;

			for ( auto iter = begin; iter != end; iter++ ) {
				if ( best < *iter ) best = *iter;
			}

			return best;
		}

		template<typename T, typename Iter, typename Func>
		static T Max( Iter begin, Iter end, Func f, T initial ) {
			T best = initial;

			for ( auto iter = begin; iter != end; iter++ ) {
				T v = f( *iter );
				if ( best < v ) best = v;
			}

			return best;
		}

		template<typename T, typename Iter>
		static T Min( Iter begin, Iter end, T initial ) {
			T best = initial;

			for ( auto iter = begin; iter != end; iter++ ) {
				if ( *iter < best ) best = *iter;
			}

			return best;
		}

		template<typename T, typename Iter, typename Func>
		static T Min( Iter begin, Iter end, Func f, T initial ) {
			T best = initial;

			for ( auto iter = begin; iter != end; iter++ ) {
				T v = f( *iter );
				if ( v < best ) best = v;
			}

			return best;
		}

		template<typename T, typename Iter>
		static T Max( Iter begin, size_t n, T initial ) {
			T best = initial;
			auto iter = begin;

			for ( size_t i = 0; i < n; i++, iter++ ) {
				if ( best < *iter ) best = *iter;
			}

			return best;
		}

		template<typename T, typename Iter, typename Func>
		static T Max( Iter begin, size_t n, Func f, T initial ) {
			T best = initial;
			auto iter = begin;

			for ( size_t i = 0; i < n; i++, iter++ ) {
				T v = f( *iter );
				if ( best < v ) best = v;
			}

			return best;
		}

		template<typename T, typename Iter>
		static T Min( Iter begin, size_t n, T initial ) {
			T best = initial;
			auto iter = begin;

			for ( size_t i = 0; i < n; i++, iter++ ) {
				if ( *iter < best ) best = *iter;
			}

			return best;
		}

		template<typename T, typename Iter, typename Func>
		static T Min( Iter begin, size_t n, Func f, T initial ) {
			T best = initial;
			auto iter = begin;

			for ( size_t i = 0; i < n; i++, iter++ ) {
				T v = f( *iter );
				if ( v < best ) best = v;
			}

			return best;
		}

		template<typename T, typename Iter>
		static T Accumulate( Iter begin, size_t n, T initial ) {
			T sum = initial;
			auto iter = begin;

			for ( size_t i = 0; i < n; i++, iter++ ) {
				sum += *iter;
			}

			return sum;
		}

		template<typename T, typename Iter, typename AddFunc>
		static T Accumulate( Iter begin, size_t n, AddFunc f, T initial ) {
			T sum = initial;
			auto iter = begin;

			for ( size_t i = 0; i < n; i++, iter++ ) {
				sum = f( *iter, sum );
			}

			return sum;
		}

		template<typename CollectionType>
		static void Free(CollectionType& collection) {
			for ( auto &p : collection ) {
				delete p;
			}
		}

		template<typename CollectionType, typename ElementType>
		static void Fill( CollectionType& collection, ElementType value ) {
			std::fill( collection.begin(), collection.end(), value );
		}

		template<typename X, typename Y>
		static void LinFit( const X& x, const Y& y, size_t n, double& a, double& b ) {
			double sumX = 0, sumY = 0, sumXY = 0, sumXSquared = 0;

			for ( size_t i = 0; i < n; i++ ) {
				if ( !std::isfinite( y[i] ) ) continue;
				sumX += x[i];
				sumY += y[i];
				sumXY += x[i] * y[i];
				sumXSquared += x[i] * x[i];
			}

			double A = sumXSquared, B = sumX, C = sumX, D = n;
			double det = A * D - B * C;
			double A1 = D / det, B1 = -B / det, C1 = -C / det, D1 = A / det;
			a = A1 * sumXY + B1 * sumY;
			b = C1 * sumXY + D1 * sumY;
		}

		template<typename X, typename Y>
		static Y Lerp( const X& x, double& a, double& b ) {
			return a * x + b;
		}

		template<typename X, typename Y>
		static Y Lerp( const X& x, const X& x0, const Y& y0, const X& x1, const Y& y1 ) {
			return x1 == x0 ? (Y) NAN : y0 + ( x - x0 ) * ( y1 - y0 ) / ( x1 - x0 );
		}

		static double LogOnePlusX( double x ) {
			if ( fabs( x ) >= 1e-10 ) {
				return log( 1 + x );
			}

			double sum = 0;
			double currentTerm = 1;
			int sign = -1;

			for ( int i = 1;; i++ ) {
				sign = -sign;
				currentTerm *= x;

				double prev = sum;

				sum += sign * currentTerm / i;

				if ( sum == prev ) {
					break;
				}

				//std::cerr << "i = " << i << ", currentTerm = " << currentTerm << ", sum = " << sum << "\n";
			}

			//std::cerr << "\n";

			return sum;
		}

		static double OneMinusExpX( double x ) {
			if ( fabs( x ) >= 1e-10 ) {
				return 1 - exp( x );
			}

			double sum = 0;
			double currentTerm = 1;

			for ( int i = 1; ; i++ ) {
				currentTerm *= x / i;

				double prev = sum;

				sum += currentTerm;

				if ( sum == prev ) break;

				//std::cerr << "i = " << i << ", currentTerm = " << currentTerm << ", sum = " << sum << "\n";
			}

			//std::cerr << "\n";

			return -sum;
		}

		static string PrintableChars() {
			string s;
			for ( uint8_t ch = 32; ch < 128; ch++ ) {
				s.push_back( ch );
			}
			return s;
		}

		template<typename T>
		static T Parse( const string& text ) {
			T val;
			istringstream str( text );
			str >> val;
			return val;
		}

		template<typename T>
		static string ToString( const T& value ) {
			ostringstream str;
			str << value;
			return str.str();
		}
	};

	class Int {
	public:
		static int Parse( const string& s ) {
			try {
				size_t read = 0;
				return std::stoi( s, &read );
			}
			catch ( std::invalid_argument& ) {
				throw FormatException( "Invalid integer data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( int value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<int>& x, string& delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Uint {
	public:
		static uint Parse( const string& s ) {
			try {
				size_t read = 0;
				return (uint) std::stoul( s, &read );
			}
			catch ( std::invalid_argument& ) {
				throw Exception( "Invalid unsigned integer data in string '" + s + "'", FileAndLine );
			}

		}

		static string ToString( uint value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<uint>& x, string& delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Uint64 {
	public:
		static uint64_t Parse( const string& s ) {
			try {
				size_t read = 0;
				return std::stoull( s, &read );
			}
			catch ( std::invalid_argument& ) {
				throw Exception( "Invalid floating point data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( uint64_t value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<uint64_t>& x, string& delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Int64 {
	public:
		static int64_t Parse( const string& s ) {
			try {
				size_t read = 0;
				return std::stoll( s, &read );
			}
			catch ( std::invalid_argument& ) {
				throw Exception( "Invalid floating point data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( int64_t value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<int64_t>& x, string& delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	template <typename T>
	class Convert {
	public:
		static T Parse( const string& s ) {
			istringstream str( s );
			T result;
			str >> result;
			return result;
		}

		static string ToString( T value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<T>& x, string& delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Ulong {
	public:
		static ulong Parse( const string& s ) {
			try {
				size_t read = 0;
				return std::stoul( s, &read );
			}
			catch ( std::invalid_argument& ) {
				throw Exception( "Invalid unsigned long data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( ulong value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<ulong>& x, string& delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Double {
		/**
		*	<summary>
		*	Scans a floating point value from a supplied string.
		*	</summary>
		*	<param name="s">The string which is isupposed to contain text
		*	representation of a floating point value.</param>
		*	<returns>
		*	</returns>
		*/
	public:
		static double Parse( const string& s ) {
			try {
				size_t read = 0;
				return std::stod( s, &read );
			}
			catch ( std::invalid_argument& ) {
				throw Exception( "Invalid floating point data in string '" + s + "'", FileAndLine );
			}
		}

		static string ToString( double value ) {
			ostringstream str;
			str << value;
			return str.str();
		}

		static string Join( const vector<double>& x, string& delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << *current;
				}
			}

			return s.str();
		}
	};

	class Bool {
	public:
		static bool Parse( const string& s ) {
			const char
				* t = "true",
				* cs = s.c_str();

			for ( ; *cs && *t; cs++, t++ ) {
				if ( tolower( *cs ) != *t ) return false;
			}

			return *cs == 0 && *t == 0;
		}

		static string ToString( double value ) {
			ostringstream str;
			str << ( value ? "true" : "false" );
			return str.str();
		}

		static string Join( const vector<bool>& x, string& delimiter ) {
			ostringstream s;

			if ( x.size() > 0 ) {
				s << x[0];

				for ( auto current = x.begin() + 1; current != x.end(); current++ ) {
					s << delimiter << ( *current ? 1 : 0 );
				}
			}

			return s.str();
		}
	};

	class File {
	public:
		static bool Exists( const string& fileName ) {
			FILE* f = fopen( fileName.c_str(), "r" );
			if ( f ) {
				fclose( f );
				return true;
			}
			return false;
		}

		static void ReadStrings( std::istream& str, std::function<void( const string& s )> action ) {
			string s;
			while ( std::getline( str, s ) ) {
				action( s );
			}
		}

		static void ReadStrings( const string& fileName, std::function<void( const string& s )> action ) {
			ifstream str( fileName );
			ReadStrings( str, action );
		}
	};
}
