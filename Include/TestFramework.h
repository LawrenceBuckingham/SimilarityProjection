#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "Delegates.hpp"
#include "Exception.hpp"

using namespace std;

namespace QutBio {

	struct TestRecord {
		string name;
		Action test;
		string desc;
	};

	class TestFramework {
		/// <summary> Executes all functions that have a name starting with the 
		/// prefix "test", which corresponds to C# names starting with "Test".
		/// <para>
		///		Each test will be executed against a freshly constructed instance 
		///		of the type.
		/// </para>
		/// </summary>
		/// <param name="types"></param>

	public:
		static void RunAllTests( vector<TestRecord> & tests ) {
			int passed = 0;
			int outOf = 0;

			for ( TestRecord test : tests ) {
				outOf++;
				try {
					test.test();
					cout << test.name << "(" << test.desc << ") passed!" << endl;
					passed++;
				}
				catch ( Exception & ex ) {
					cout << test.name << "(" << test.desc << ") failed!" << endl;
					cout << ex.File() << ":" << ex.Line() << " - " << ex.what() << endl;
				}
			}

			cout << "Passed " << passed << "/" << outOf << ", failed " << ( outOf - passed ) << "/" << outOf << endl << endl;
		}

	};
}
