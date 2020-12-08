#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <iostream>
#include <string>
#include <unordered_map>

#include "Exception.hpp"
#include "CsvIO.hpp"

namespace QutBio{
	class FileUtil{
	public:
		static void WriteCsvRecord( ostream & outFD, vector<string> record ){
			CsvWriter io( outFD );
			io.Write( record );
		}
	};
}
