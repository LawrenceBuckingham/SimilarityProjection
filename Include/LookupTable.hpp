#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <unordered_map>

#include "Assert.hpp"
#include "Delegates.hpp"

using namespace QutBio;

namespace QutBio {
	/**
	<summary>
		A lookup table mapping some key of type K to _pointer_to_ values of type V, which owns the stored pointers.
		The values are dynamically allocated when a value is added by a factory function which is responsible for
		correct initialisation.
		When the collection is deleted the pointers are destroyed.
	<summary>
		*/

	template<typename K, typename V>
	class LookupTable {
	private:
		unordered_map<K, V *> map;
	public:
		LookupTable() {}

		~LookupTable() {
			for ( auto kv : map ) {
				delete kv.second;
			}
		}

		LookupTable( const LookupTable & other ) = delete;

		LookupTable & operator = ( const LookupTable & other ) = delete;

		void Add( const K & key, Func<V *> factory ) {
			assert_false( map.find( key ) != map.end() );
			V * value = factory();
			map[key] = value;
		}

		V * operator[] ( const K & key ) {
			return map[key];
		}

		void ForEach( Action2<K&, V *> action ) {
			for ( auto kv : map ) {
				action( kv.first, kv.second );
			}
		}
	};
}
