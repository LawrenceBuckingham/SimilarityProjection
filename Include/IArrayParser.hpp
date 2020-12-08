#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>
#include <string>


namespace QutBio{

	class IArrayParser {
	public: virtual void Parse( vector<string> & fields ) = 0;
	};

}
