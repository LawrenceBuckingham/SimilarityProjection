#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <iostream>

namespace QutBio{
	struct TrecEvalRecord {
		string queryId;
		string subjectId;
		double similarity;

		friend ostream & operator << (ostream & out, const TrecEvalRecord & trec) {
			return (out << trec.queryId << " 0 " << trec.subjectId << " 0 " << trec.similarity << " ignored");
		}
	};
}
