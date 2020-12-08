#pragma once

#include <vector>

namespace QutBio {
	template<typename T>
	class Pointers {
	private:
		std::vector<T *> items;

	public:
		virtual ~Pointers() { for ( auto i : items ) delete i; }

		T * operator()() {
			items.push_back( new T() );
			return items.back();
		}

		template<typename ...Args>
		T *operator()( Args... args ) {
			T *t =  new T( args... );
			items.push_back(t);
			return t;
		}

	};

}
