#pragma once
#include <unordered_map>
#include <vector>

namespace QutBio {
	template<typename KeyType, typename HashType>
	class Registry {
		unordered_map<KeyType, size_t, HashType> reg;
		vector<KeyType> key;

	public:
		/// <summary>
		/// Gets a canonical id number corresponding to a key.
		/// </summary>
		/// <param name="id"></param>
		/// <returns></returns>
		size_t operator()( const KeyType & id ) {
			auto i = reg.find( id );

			if ( i == reg.end() ) {
				size_t n = reg.size();
				pair<KeyType, size_t> p{ id, n };
				reg.insert( p );
				key.push_back(id);
				return n;
			}
			else {
				return i->second;
			}
		}

		size_t Size() const {
			return reg.size();
		}

		auto Begin() const {
			return reg.begin();
		}

		auto End() const {
			return reg.end();
		}

		void Clear() {
			reg.clear();
		}

		const KeyType& At( size_t i ) {
			return key.at(i);
		}
	};
}
