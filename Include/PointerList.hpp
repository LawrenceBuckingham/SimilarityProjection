#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>

#include "Util.hpp"

namespace QutBio {
	template<typename T>
	class PointerList : public std::vector<T *> {
		using BaseClass = std::vector<T *>;
	public:
		PointerList(const PointerList & other) = delete;

		PointerList(size_t size = 0) : BaseClass(size) {}

		virtual ~PointerList() {
			Util::Free(*this);
		}
	};
}
