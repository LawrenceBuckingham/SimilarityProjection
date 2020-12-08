#pragma once

// Trick Visual Studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <algorithm>
#include <vector>
#include <iostream>
#include <string.h>

using namespace std;

#include "db.hpp"
#include "Delegates.hpp"
#include "Exception.hpp"
#include "Types.hpp"

namespace QutBio {

	/**
	 *	The idea here is to get multiple views on the same underlying array.
	 *	I will experiment with an object-based approach (Version 1) and
	 *	when I'm sure I have the logic correct convert it to something more
	 *	efficient (hopefully) (pack a list of signatures into a contiguous extent of
	 *	memory to improve cache locality).
	 */

	template<typename ItemType>
	class PackedArray {
	protected:
		typedef uint64_t Word;
		const uint digits = 64;
		size_t bitsPerItem;
		size_t length;
		vector<Word> items;

	public:
		PackedArray(
			size_t bitsPerItem,
			size_t length
		) :
			bitsPerItem( bitsPerItem ),
			length( length )
			//
		{
			if ( bitsPerItem > digits ) {
				throw Exception( "Invalid bitsPerItem for packed array.", FileAndLine );
			}

			// Get the number of words required, plus one as a zero-filled
			// overflow region at the end.
			size_t requiredCapacity = ( length * bitsPerItem + digits - 1 ) / digits;
			items.resize( requiredCapacity + 1, 0ull );
		}

		PackedArray( 
			const string & s, 
			size_t bitsPerSymbol, 
			const vector<uint8_t> & inverse
		) : PackedArray( bitsPerSymbol, s.length() ) {
			for ( size_t i = 0; i < s.length(); i++ ) {
				int item = inverse[(uint8_t) s[i]];
				Set( i, item );
			}
		}

		template<typename Container>
		PackedArray(
			const Container & s,
			size_t bitsPerSymbol,
			const byte inverse[]
		) : PackedArray(bitsPerSymbol, s.size()) {
			for (size_t i = 0; i < s.size(); i++) {
				int item = inverse[(uint8_t)s[i]];
				Set(i, item);
			}
		}

		PackedArray(
			const char * s,
			size_t size,
			size_t bitsPerSymbol,
			const byte inverse[]
		) : PackedArray(bitsPerSymbol, size) {
			for (size_t i = 0; i < size; i++) {
				int item = inverse[(uint8_t)s[i]];
				Set(i, item);
			}
		}

		template<typename Container>
		PackedArray( 
			const Container & s, 
			size_t bitsPerSymbol 
		) : PackedArray( bitsPerSymbol, s.size() ) {
			for ( size_t i = 0; i < s.size(); i++ ) {
				ItemType item = s[i];
				Set( i, item );
			}
		}

		PackedArray( const PackedArray & other ) :
			bitsPerItem( other.bitsPerItem ),
			length( other.length ),
			items( other.items )
			//
		{}

		PackedArray & operator=( const PackedArray & other ) {
			bitsPerItem = other.bitsPerItem;
			length = other.length;
			items = other.items;
		}

		virtual ~PackedArray() {}

		size_t Length() const {
			return length;
		}

		size_t BitsPerItem() const {
			return bitsPerItem;
		}

		ItemType Get( size_t index ) const {
			uint pos = index * bitsPerItem;
			uint wordIndex = pos / digits;
			uint offset = pos % digits;
			size_t n = items.size() - 1;

			if ( wordIndex >= n ) throw Exception( "Index out of bounds", FileAndLine );

			if ( offset <= digits - bitsPerItem ) {
				Word mask = ( Word( 1 ) << bitsPerItem ) - 1;
				return ( items[wordIndex] >> offset ) & mask;
			}
			else {
				uint overlap = digits - offset;
				Word mask = ( ( Word( 1 ) << bitsPerItem ) - 1 ) >> overlap;
				return ( items[wordIndex] >> offset ) | ( ( items[wordIndex + 1] & mask ) << overlap );
			}
		}

		PackedArray & Set( size_t index, ItemType value ) {
			uint pos = index * bitsPerItem;
			uint wordIndex = pos / digits;
			uint offset = pos % digits;
			size_t n = items.size() - 1;

			if ( wordIndex >= n ) throw Exception( "Index out of bounds", FileAndLine );

			if ( offset <= digits - bitsPerItem ) {
				items[wordIndex] |= Word( value ) << offset;
			}
			else {
				uint overlap = digits - offset;
				Word mask = ( Word( 1 ) << overlap ) - 1;
				items[wordIndex] |= ( mask & value ) << offset;
				items[wordIndex + 1] = value >> overlap;
			}

			return *this;
		}

		bool operator==( const PackedArray & other ) const {
			if ( other.length != length ) {
				throw Exception( "length does not match", FileAndLine );
			}

			const uint n = items.size();

			for ( int i = 0; i < n; i++ ) {
				if ( items[i] != other.items[i] ) {
					return false;
				}
			}

			return true;
		}

		bool operator!=( const BitSet & other ) const {
			return !operator==( other );
		}

		template<typename T>
		void Unpack( vector<T> & dest ) {
			for ( uint i = 0; i < length; i++ ) {
				dest.push_back( Get( i ) );
			}
		}

		void Unpack( string & dest, const string & symbols ) {
			for ( uint i = 0; i < length; i++ ) {
				dest.push_back( symbols[Get( i )] );
			}
		}

		size_t BitsPerItem() {
			return bitsPerItem;
		}

		size_t Length() {
			return length;
		}

		const vector<Word> & Items() {
			return items;
		}

		friend ostream & operator << (ostream & s, const PackedArray & x) {
			size_t len = 0;
			size_t max = x.length * x.bitsPerItem;

			for ( auto word: x.items ) {
				for (size_t i = 0; len < max && i < 64; i++, len++) {
					s << ((word >> i) & 1);
				}
			}

			return s;
		}

		Word GetSlice(size_t slice, size_t slices) {
			size_t totalBits = length * bitsPerItem;
			size_t firstBit = slice * totalBits / slices;
			size_t lastBit = (slice + 1) * totalBits / slices - 1;

#define bitsPerItem  ERROR!!!
			size_t bitsPerSlice = lastBit - firstBit + 1;

			if ( bitsPerSlice > digits ) {
				throw Exception( "Slice width exceeds size of word", FileAndLine );
			}

			size_t firstWord = firstBit / digits;
			size_t firstOffset = firstBit % digits;

			if (firstOffset <= digits - bitsPerSlice) {
				Word mask = (Word(1) << bitsPerSlice) - 1;
				return (items[firstWord] >> firstOffset) & mask;
			}
			else {
				uint overlap = digits - firstOffset;
				Word mask = ((Word(1) << bitsPerSlice) - 1) >> overlap;
				return (items[firstWord] >> firstOffset) | ((items[firstWord + 1] & mask) << overlap);
			}
#undef bitsPerItem
		}
	};

}
