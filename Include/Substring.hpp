#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cstddef>
#include <cstring>
#include <cstdio>
#include <ios>
#include <functional>

namespace QutBio {
	// Class representing a sub-string of a raw char array. 
	class Substring {
	private:
		// The address of the first character in the instance. 
		// NB: strlen(string) > length with high probability
		const Symbol* chars;

		// The number of characters in the kmer.
		size_t length;

		// Cached hash code.
		size_t hashCode;

	public:
		Substring(
			const Symbol* str,
			size_t start, 
			size_t length,
			const Alphabet * alphabet = 0
			) : 
				chars(str + start), 
				length(length), 
				hashCode( alphabet ? alphabet->Encode<size_t>(chars, length) : HashCode( chars, length ) ) {}

		Substring(const Substring & other) {
			this->chars = other.chars;
			this->length = other.length;
			this->hashCode = other.hashCode;
		}

		Substring & operator=(const Substring & other) {
			this->chars = other.chars;
			this->length = other.length;
			this->hashCode = other.hashCode;
			return *this;
		}

		virtual ~Substring() {}

		const Symbol* Chars() const { return chars; }

		size_t Length() const { return length; }

		size_t size() const { return length; }

		const Symbol& operator[](size_t index) const {
			return chars[index];
		}

		friend bool operator==(const Substring & lhs, const Substring & rhs) {
			return memcmp(lhs.chars, rhs.chars, lhs.length) == 0;
		}

		friend bool operator!=(const Substring & lhs, const Substring & rhs) {
			return memcmp(lhs.chars, rhs.chars, lhs.length) != 0;
		}

		friend bool operator<(const Substring & lhs, const Substring & rhs) {
			return memcmp(lhs.chars, rhs.chars, lhs.length) < 0;
		}

		friend std::ostream & operator<<(std::ostream & str, const Substring & s) {
			for (size_t i = 0; i < s.length; i++) {
				str << s.chars[i];
			}

			return str;
		}

		struct Hash {
			size_t operator()(const Substring & __val) const noexcept {
				return __val.HashCode();
			}
		};

		size_t HashCode() const {
			return hashCode;
		}

		static size_t HashCode( const Symbol* chars, size_t length ) {
			return _Hash_impl::hash( (void *) chars, length );
		}

		void fprint( FILE * stream ) const {
			for (size_t i = 0; i < length; i++) {
				fputc( chars[i].value, stream );
			}
		}
	};

	typedef Substring * pSubstring;
}
