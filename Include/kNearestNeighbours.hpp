#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>
using namespace std;

namespace QutBio {
	/**
	*  <p>
	*      Uses a fixed-size min heap to implement a k-nearest neighbour
	*      accumulator.
	*  </p>
	*  <p>
	*      The data is stored in approximately ascending order (the heap
	*      property) with the greatest value at the back of the list. (Usually,
	*      a min heap would have the smallest value at the front but the
	*      greatest value could be some place near, but not necessarily
	*      at, the back.)
	*       The object has a fixed maximum capacity K, and is organised so that
	*       the K nearest neighbours are stored (which is why the largest value
	*       must be at the back). The items are extracted and sorted into true
	*       ascending order afterwards.
	*  </p>
	*/
	template<class T, typename Compare = less<T>>
	class KnnHeap {
	protected:
		vector<T> heap;
		const Compare & compare;
		size_t capacity;
	public:
		KnnHeap(size_t capacity = 0, const Compare & compare = Compare()) : compare(compare), capacity(capacity) {
			heap.reserve(capacity);
		}

		void setCapacity(size_t capacity) {
			this->capacity = capacity;
			heap.reserve(capacity);
		}

		size_t getCapacity() {
			return capacity;
		}

		void clear() {
			heap.clear();
		}

		void push(const T & item) {
			size_t last = heap.size() - 1;

			if (last + 1 == capacity) {
				if (compare(item, heap[last])) {
					heap[last] = item;
				}
			}
			else {
				heap.push_back(item);
			}

			make_heap(heap.rbegin(), heap.rend(), compare);
		}

		template<typename Collection>
		void push(const Collection & kmers) {
			for (T item : kmers) {
				push(item);
			}
		}

		const T & top() const {
			return heap.back();
		}

		void pop() {
			heap.pop_back();
			make_heap(heap.rbegin(), heap.rend(), compare);
		}

		bool empty() {
			return heap.empty();
		}

		typedef typename std::vector<T>::iterator iterator;

		iterator begin() {
			auto iter = heap.begin();
			return iter;
		}

		iterator end() {
			auto iter = heap.end();
			return iter;
		}

		friend ostream & operator << (ostream & out, KnnHeap<T, Compare> & knn) {
			bool deja = false;

			for (auto x : knn) {
				if (deja) cout << ", "; else deja = true;
				out << x;
			}

			out << endl;
			return out;
		}
	};

	/**
	*  <p>
	*      Uses a fixed-size unordered vector to implement a k-nearest neighbours list.
	*  </p>
	*  <p>
	*       The object has a fixed maximum capacity K, and is organised so that
	*       the K "smallest" items are stored.
	*		The data is stored in haphazard order - initially the order in which the first
	*		K items arrive. The location of the largest item in the collection is saved.
	*		Once the capacity K is reached, the largest item removed if a smaller item is
	*		inserted, and the container is scanned once again to update the location of
	*		the largest.
	*		The collection can be sorted into true ascending order at any time by calling
	*		Sort.
	*  </p>
	*/
	template<typename ElementType, typename Distance>
	struct KnnVector {
		using Record = pair<Distance, ElementType>;
		vector<Record> elements;
		size_t capacity;
		Distance minDistance;
		Distance ejectDistance = numeric_limits<Distance>::min();
		size_t ejectPos = 0;

		KnnVector(size_t capacity, Distance minDistance) : capacity(capacity), minDistance(minDistance), ejectDistance(minDistance) {
			elements.reserve(capacity);
		}

		void setCapacity(size_t capacity) {
			this->capacity = capacity;
			elements.reserve(capacity);
		}

		size_t getCapacity() {
			return capacity;
		}

		void clear() {
			elements.clear();
			ejectPos = 0;
			ejectDistance = minDistance;
		}

		bool canPush(const Distance & distance) {
			return elements.size() < capacity || distance < ejectDistance;
		}

		void push(const ElementType & item, const Distance & distance) {
			if (elements.size() < capacity) {
				if (distance > ejectDistance) {
					ejectDistance = distance;
					ejectPos = elements.size();
				}

				elements.emplace_back(distance, item);
			}
			else if (distance < ejectDistance) {
				elements[ejectPos].first = distance;
				elements[ejectPos].second = item;

				ejectPos = 0;
				ejectDistance = elements[0].first;

				for (uint i = 1; i < capacity; i++) {
					Distance d = elements[i].first;

					if (d > ejectDistance) {
						ejectPos = i;
						ejectDistance = d;
					}
				}
			}
		}

		bool empty() {
			return elements.empty();
		}

		void sort() {
			std::sort( elements.begin(), elements.end(), 
				[] ( const Record & lhs, const Record & rhs ) {
					return lhs.first < rhs.first;
				}
			);
		}

		typedef typename std::vector<Record>::iterator iterator;

		iterator begin() {
			auto iter = elements.begin();
			return iter;
		}

		iterator end() {
			auto iter = elements.end();
			return iter;
		}

		friend ostream & operator << (ostream & out, KnnVector<ElementType, Distance> & knn) {
			bool deja = false;

			for (auto & x : knn) {
				if (deja) cout << ", "; else deja = true;
				out << x.second;
			}

			out << endl;
			return out;
		}
	};
}
