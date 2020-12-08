#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <vector>

#include "FastaSequence.hpp"

namespace QutBio {
	struct SequenceWrapper {
		FastaSequence * base;

		SequenceWrapper(FastaSequence * base) : base(base) {}

		const string & IdStr() {
			return base->IdStr();
		}

		size_t KmerCount(uint k) {
			return base->KmerCount(k);
		}

		const vector<Symbol> & Sequence() {
			return base->Sequence();
		}

		const string & CharData() {
			return base->CharData();
		}

		size_t Length() const {
			return base->Length();
		}

		string DefLine() const {
			return base->DefLine();
		}

		template<typename T, typename C>
		static void Wrap( const C & collection, vector<T *> result) {
			for ( auto seq: collection ) {
				result.push_back( new T( seq ) );
			}
		}

		friend ostream &operator<<(ostream &out, const SequenceWrapper * seq) {
			out << seq->base;
			return out;
		}

	};
}
