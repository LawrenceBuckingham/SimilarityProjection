#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <iostream>
#include <string>

namespace QutBio {
	class Console {
	public:
		static void Error( const string & s ) {
			cerr << s << endl;
		}

		static void Log( const string & s ) {
			cout << s << endl;
		}

		static void Error( const char * s ) {
			cerr << s << endl;
		}

		static void Log( const char s ) {
			cout << s << endl;
		}
	};
}
