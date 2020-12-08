#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <algorithm>
#include <vector>
#include <iostream>
#include <string.h>

//#if USE_OMP
//#include <omp.h>
//#define LOCK omp_set_lock(&lock)
//#define UNLOCK omp_unset_lock(&lock)
//#define // LOCK_WORD(expr) omp_set_lock(&locks[expr])
//#define // UNLOCK_WORD(expr) omp_unset_lock(&locks[expr])
//#else
//#define LOCK
//#define UNLOCK
//#define // LOCK_WORD(expr)
//#define // UNLOCK_WORD(expr)
//#endif

using namespace std;

#include "db.hpp"
#include "Delegates.hpp"
#include "Exception.hpp"

namespace QutBio
{

/*
	**	<summary>
	**		Set of 0..(capacity-1), implemented via bit-packed array.
	**	</summary>
	*/

class BitSet: protected vector<uint64_t>
{
  protected:
	typedef uint64_t Bits;
	static const size_t DIGITS = 64;
	static const size_t SHIFT = 6;
	static const size_t MASK = (1 << SHIFT) - 1;
	size_t capacity_;
	//#if USE_OMP
	//		vector<omp_lock_t> locks;
	//		omp_lock_t lock;
	//#endif

  public:
	/*
		**	Summary:
		**		Construct a BitSet.
		**		Complexity: O(capacity).
		**	Parameters:
		**		capacity - the number of bits in the set. Valid members are therefore
		**				0..(capacity-1).
		*/
	BitSet(size_t capacity = DIGITS) : vector<uint64_t>((capacity + DIGITS - 1) >> SHIFT),
									   capacity_(capacity)
	//#if USE_OMP
	//			,
	//			locks((capacity + DIGITS - 1) >> SHIFT)
	//#endif
	//
	{
		std::fill(begin(), end(), 0);
	}

	/*
		**	Summary:
		**		Construct a BitSet by copying another.
		**		Complexity: O(capacity).
		**	Parameters:
		**		capacity - the number of bits in the set. Valid members are therefore
		**				0..(capacity-1).
		*/
	BitSet(const BitSet &other) : vector<uint64_t>(other), capacity_(other.capacity_)
	//
	{
		//#if USE_OMP
		//			for (auto & lock : locks) {
		//				omp_init_lock(&lock);
		//			}
		//
		//			omp_init_lock(&lock);
		//#endif
	}

	/*
		**	Summary:
		**		Copy a BitSet from another.
		**		Complexity: O(capacity).
		**
		**	Decided that I don't want this.
		**
		**	Parameters:
		**		capacity - the number of bits in the set. Valid members are therefore
		**				0..(capacity-1).
		*/
	BitSet &operator=(const BitSet &other)
	{
		// LOCK;
		vector<uint64_t>::operator=(other);
		capacity_ = other.capacity_;
		return *this;
		//UNLOCK;
	}

	/*
		**	Summary:
		**		Destroy a BitSet.
		*/
	virtual ~BitSet()
	{
		//#if USE_OMP
		//			omp_destroy_lock(&lock);
		//
		//			for (auto &l : locks) {
		//				omp_destroy_lock(&l);
		//			}
		//#endif
	}

	/*
		**	Summary:
		**		Get the number of elements that can be stored in the set.
		**		Complexity: O(1).
		*/
	size_t Capacity() const
	{
		return capacity_;
	}

	/*
		**	Summary:
		**		Get the number of elements that are currently stored in the set.
		**		Complexity: O(C), where C is capacity.
		*/
	size_t Cardinality()
	{
		// LOCK;

		size_t cardinality = 0;

		for (auto w : *this)
		{
			cardinality += POPCOUNT(w);
		}

		//UNLOCK;

		return cardinality;
	}

	/*
		Summary:
			Determine whether set contains a designated value.
			Complexity: O(1)
		Parameters:
			i - the value to locate.
		Returns:
			true iff i is in the current set.
		*/
	bool Contains(size_t i)
	{
		if (i >= capacity_)
		{
			throw Exception("Index out of bounds", FileAndLine);
		}
		// LOCK_WORD(i >> SHIFT);
		bool result = (*this)[i >> SHIFT] & (((Bits)1) << (i & MASK));
		// UNLOCK_WORD(i >> SHIFT);
		return result;
	}

	/*
		Summary:
			Determine whether set is empty.
			Complexity: O(C), where C is capacity.
		Returns:
			true iff cardinality of the current set is zero.
		*/
	bool IsEmpty()
	{
		return Cardinality() == 0;
	}

	/*
		**	Summary:
		**		Insert a value into the set.
		**		Complexity: O(1).
		**	Parameters:
		**		i - a value in 0..(capacity-1).
		*/
	void Insert(size_t i)
	{
		if (i >= capacity_)
		{
			throw Exception("Index out of bounds", FileAndLine);
		}

		// LOCK_WORD(i >> SHIFT);
		(*this)[i >> SHIFT] |= (((Bits)1) << (i & MASK));
		// UNLOCK_WORD(i >> SHIFT);
	}

	/*
		**	Summary:
		**		Insert a set of values into the current set.
		**		Complexity: O(s.Capacity())
		**	Parameters:
		**		s - A set having capacity equal to or less than
		**			the current set.
		*/
	void Union(const BitSet &s)
	{

		if (s.capacity_ > capacity_)
		{
			throw Exception("Inserted Bit Set is too large", FileAndLine);
		}

		// LOCK;

		for (size_t i = 0; i < size(); i++)
		{
			(*this)[i] |= s[i];
		}

		//UNLOCK;
	}

	/*
		**	Warning:
		**		Uses omp parallel for, so don't call this from within a parallel block!
		**	Summary:
		**		Insert a set of values into the current set.
		**		Complexity: O(s.Capacity())
		**	Parameters:
		**		s - A set having capacity equal to or less than
		**			the current set.
		*/
	void Union(const vector<BitSet> &s)
	{
		// LOCK;
		for (size_t i = 0; i < size(); i++)
		{
			for (size_t j = 0; j < s.size(); i++)
			{
				(*this)[i] |= s[j][i];
			}
		}
		//UNLOCK;
	}

	/*
		**	Summary:
		**		Insert a list of values into the current set
		**		Complexity: O(s.size())
		**	Parameters:
		**		i - a list of values in 0 .. (capacity-1).
		*/
	template <typename T>
	void Union(const vector<T> &s)
	{
		// LOCK;
		for (auto i : s)
		{
			Insert(i);
		}
		//UNLOCK;
	}

	/*
		**	Summary:
		**		Remove a value from the set. If the designated value is
		**		not present in the set, then nothing happens.
		**		Complexity: O(1)
		**	Parameters:
		**		i - a value in 0..(capacity-1).
		*/
	void Remove(size_t i)
	{
		if (i >= capacity_)
		{
			throw Exception("Index out of bounds", FileAndLine);
		}

		// UNLOCK_WORD(i >> SHIFT);
		(*this)[i >> SHIFT] &= ~(((Bits)1) << (i & MASK));
		// LOCK_WORD(i >> SHIFT);
	}

	/*
		**	Summary:
		**		Remove a set (s) of values from the current set.
		**		Complexity: O(s.capacity).
		**	Parameters:
		**		s - a set having capacity equal to or less than that
		**			of the current set.
		*/
	void Remove(const BitSet &s)
	{
		// LOCK;
		for (size_t i = 0; i < s.size(); i++)
		{
			(*this)[i] &= ~(s[i]);
		}
		//UNLOCK;
	}

	/*
		**	Summary:
		**		Insert a list of values into the current set
		**		Complexity: O(s.size())
		**	Parameters:
		**		i - a list of values in 0 .. (capacity-1).
		*/
	template <typename T>
	void Remove(const vector<T> &s)
	{
		// LOCK;
		for (auto i : s)
		{
			Remove(i);
		}
		//UNLOCK;
	}

	/*
		**	Summary:
		**		Flips all bits in the current set.
		**		Complexity: O(capacity)
		**	Returns:
		**		A reference to the current set.
		*/
	BitSet &Complement()
	{
		// LOCK;

		for (auto &word : *this)
		{
			word = ~word;
		}

		uint64_t bitsInLastWord = capacity_ & 63ull;

		if (bitsInLastWord)
		{
			this->back() &= (1ull << bitsInLastWord) - 1ull;
		}

		//UNLOCK;

		return *this;
	}

	/*
		**	Summary:
		**		Determine if the current set is equal to another set.
		**		Complexity: O(capacity)
		**	Returns:
		**		A reference to the current set.
		*/
	bool operator==(const BitSet &other)
	{
		bool result = true;

		if (other.capacity_ != capacity_)
		{
			throw Exception("Capacity does not match", FileAndLine);
		}

		// LOCK;
		for (size_t i = 0; i < size(); i++)
		{
			if ((*this)[i] != other[i])
			{
				result = false;
			}
		}
		//UNLOCK;

		return result;
	}

	/*
		**	Summary:
		**		Determine if the current set is not equal to another set.
		**		Complexity: O(capacity)
		**	Returns:
		**		A reference to the current set.
		*/
	bool operator!=(const BitSet &other)
	{
		bool result = true;

		// LOCK;
		if (other.capacity_ != capacity_)
		{
			throw Exception("Capacity does not match", FileAndLine);
		}

		for (size_t i = 0; i < size(); i++)
		{
			if ((*this)[i] == other[i])
			{
				result = false;
			}
		}

		//UNLOCK;
		return result;
	}

	/*
		**	Summary:
		**		Intersects this set with another.
		**		Complexity: O(s.Capacity())
		**	Parameters:
		**		s - A set having capacity equal to that of
		**			the current set.
		*/
	void Intersect(const BitSet &s)
	{
		if (s.capacity_ != capacity_)
		{
			throw Exception("Inserted Bit Set is too large", FileAndLine);
		}

		// LOCK;
		for (size_t i = 0; i < s.size(); i++)
		{
			(*this)[i] &= s[i];
		}
		//UNLOCK;
	}

	void Foreach(function<void(size_t i)> callback)
	{
		size_t i = 0;

		// LOCK;
		for (auto b : *this)
		{
			if (b)
			{
				Bits mask = 1;

				for (uint j = 0; j < DIGITS; j++)
				{
					if (b & mask)
					{
						callback(i * DIGITS + j);
					}

					mask <<= 1;
				}
			}

			i++;
		}
		//UNLOCK;
	}

	//		void ParallelForeach(function<void(size_t i, uint threadId)> callback) {
	//			const size_t n = bits.size();
	//
	//#pragma omp parallel for
	//			for (size_t i = 0; i < n; i++) {
	//				uint threadId = omp_get_thread_num();
	//				uint64_t b = bits[i];
	//
	//				if (b) {
	//					Bits mask = 1;
	//
	//					for (uint j = 0; j < DIGITS; j++) {
	//						if (b & mask) {
	//							callback(i * DIGITS + j, threadId);
	//						}
	//
	//						mask <<= 1;
	//					}
	//				}
	//			}
	//		}

	void Clear()
	{
		// LOCK;
		std::fill(begin(), end(), 0ull);
		//UNLOCK;
	}

	/**
	 * Computes the similarity between this bit set and another.
	 * Similarity(x,y) = sum( i:dom(x); popcount(x[i]&y[i]) ).
	 */

	double Similarity(
		const BitSet &other) const
	{
		uint s = 0, t = 0;
		const uint n = size();

		for (uint i = 0; i < n; i++)
		{
			s += POPCOUNT((*this)[i] & other[i]);
			t += POPCOUNT((*this)[i] | other[i]);
		}

		return t == 0 ? 0.0 : double(s) / double(t);
	}

	uint HammingDistance(
		const BitSet &other)
	{
		uint s = 0;
		const uint n = size();

		for (uint i = 0; i < n; i++)
		{
			s += POPCOUNT((*this)[i] ^ other[i]);
		}

		return s;
	}

	void SetCapacity(uint capacity)
	{
		this->resize((capacity + DIGITS - 1) >> SHIFT);
		this->capacity_ = capacity;
	}

	friend ostream &operator<<(ostream &str, BitSet &bitSet)
	{
		bool deja = false;

		str << bitSet.Cardinality() << " ";

		bitSet.Foreach([&](size_t i) { str << (deja ? " " : "") << i; deja = true; });

		str << ";";

		return str;
	}

	friend istream &operator>>(istream &str, BitSet &bitSet)
	{
		bitSet.Clear();
		uint64_t t;
		size_t cardinality = 0;

		str >> cardinality;

		if (cardinality > 0)
		{
			while (!(str.eof() || str.peek() == ';'))
			{
				str >> t;
				if (str.fail())
					break;
				bitSet.Insert(t);
			}

			if (str.peek() == ';')
				str.ignore(1);
		}
		else
		{
			char c = str.get();

			while (c != ';' && !str.eof())
			{
				c = str.get();
			}
		}

		if (bitSet.Cardinality() != cardinality)
		{
			stringstream s;
			s << "bitSet cardinality " << bitSet.Cardinality() << " does not match expected value: " << cardinality;
			throw Exception(s.str(), FileAndLine);
		}

		return str;
	}
};

} // namespace QutBio
