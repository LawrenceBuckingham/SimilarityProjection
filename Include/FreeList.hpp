#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>

namespace QutBio {
	template<typename T>
	class FreeList {
		vector<T*> allItems;
		vector<T*> freeItems;
	public:
		FreeList() {}

		virtual ~FreeList() {
			for ( auto t : allItems ) {
				delete t;
			}
		}

		T*Allocate() {
			if ( freeItems.size() > 0 ) {
				auto item = freeItems.back();
				freeItems.pop_back();
				return item;
			}
			else {
				auto item = new T();
				allItems.push_back(item);
				return item;
			}
		}

		T*Allocate( function<T*(void)> factory ) {
			if ( freeItems.size() > 0 ) {
				auto item = freeItems.back();
				freeItems.pop_back();
				return item;
			}
			else {
				auto item = factory();
				allItems.push_back(item);
				return item;
			}
		}

		void Free(T* item) {
			freeItems.push_back( item );
		}

		void FreeAll() {
			freeItems = allItems;
		}
	};
}
