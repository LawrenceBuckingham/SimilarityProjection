#pragma once


// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <unordered_map>
#include <vector>
#include <string>

namespace QutBio {
	using namespace std;

	template <typename ElementType>
	class Index : public unordered_map<string, ElementType *> {
	public:
		Index () {}
		
		virtual ~Index() {}

		/**
		*	<summary>
		*		Gets an index which can be used to look up sequences by Id.
		*	</summary>
		*	<param name="db">
		*		Reference to a pointer-list of FastaSequence objects in which symbol occurrences are to be calculated.
		*	</param>
		*	<param name="index">
		*		Reference to a hash table that will be populated with the lookup table. Any previous contents in this
		*		data structure will be erased.
		*	</param>
		*/
		template <typename CollectionType>
		Index(const CollectionType & dataset) {
			for (auto seq : dataset) {
				(*this)[seq->IdStr()] = seq;
			}
		}
	};

	template <typename KeyType, typename ElementType>
	class LookupTable_ : public unordered_map<KeyType, ElementType *> {
	public:
		LookupTable_() {}
		
		virtual ~LookupTable_() {}

		/**
		*	<summary>
		*		Gets an index which can be used to look up sequences by Id.
		*	</summary>
		*	<param name="db">
		*		Reference to a pointer-list of FastaSequence objects in which symbol occurrences are to be calculated.
		*	</param>
		*	<param name="index">
		*		Reference to a hash table that will be populated with the lookup table. Any previous contents in this
		*		data structure will be erased.
		*	</param>
		*/
		template <typename CollectionType>
		LookupTable_(CollectionType & dataset) {
			for (auto seq : dataset) {
				(*this)[seq->Id()] = seq;
			}
		}
	};

}
