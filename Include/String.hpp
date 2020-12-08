#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <stddef.h>

#include <sstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cstdarg>
#include <cstdio>

#include "Types.hpp"
#include "Exception.hpp"


using namespace std;

namespace QutBio {

	class String {
	public:
		/**
		*	<summary>
		*	Creates a new copy of the supplied string in which all upper-case
		*	characters have been replaced by their lowercase equivalent, and all
		*	other characters remain unchanged.
		*	<para>
		*	The original string is unaffected.
		*	<para>
		*	<summary>
		*/
		static string ToLowerCase( const string & s ) {
			size_t len = s.length();
			string result;

			for ( size_t i = 0; i < len; i++ ) {
				auto c = s[i];

				result += tolower( c );
			}

			return result;
		}
		/**
		*	<summary>
		*	Updates the supplied string, converting all characters to their lowercase equivalent.
		*	<summary>
		*/
		static void ToLowerInPlace( string & s ) {
			std::transform( s.begin(), s.end(), s.begin(), ::tolower );
		}

		/**
		*	<summary>
		*	Creates a new copy of the supplied string in which all upper-case
		*	characters have been replaced by their lowercase equivalent, and all
		*	other characters remain unchanged.
		*	<para>
		*	The original string is unaffected.
		*	<para>
		*	<summary>
		*/
		static string ToLowerCase( const char * s ) {
			string result;

			for ( int i = 0; s[i] != 0; i++ ) {
				char c = s[i];

				result += c >= 'A' && c <= 'Z' ? c + 'a' - 'A' : c;
			}

			return result;
		}

		/** Splits a string into a sequence of tokens, breaking at the supplied delimiter character. */
		static std::vector<string>  Split( const string & s, const char * separators ) {
			std::vector<string> result;

			size_t pos = 0;

			for ( ;; ) {
				size_t index = s.find_first_of( separators, pos );

				size_t len = (index == string::npos) ? (s.length() - pos) : (index - pos);

				if ( len > 0 ) {
					result.push_back( s.substr( pos, len ) );
				}

				if ( index == string::npos ) break;

				pos = index + 1;
			}

			return result;
		}

		/** Splits a string into a sequence of tokens, breaking at the supplied delimiter character, and appending the results to a given sequence of strings. */
		static std::vector<string> Split( const string & s, char separator ) {
			std::vector<string> result;

			size_t pos = 0;

			for ( ;; ) {
				size_t index = s.find_first_of( separator, pos );

				size_t len = (index == string::npos) ? (s.length() - pos) : (index - pos);

				if ( len > 0 ) {
					result.push_back( s.substr( pos, len ) );
				}

				if ( index == string::npos ) break;

				pos = index + 1;
			}

			return result;
		}

		/** Returns a copy of a string with all characters matching the cctype isblank predicate removed.*/
		static string Trim( const string & s ) {
			string result;
			size_t len = s.size();

			if ( len > 0 ) {
				int startPos = -1;
				int endPos = -1;

				for ( size_t i = 0; i < len; i++ ) {
					if ( startPos < 0 && !isblank( s[i] ) ) {
						startPos = (int) i;
						break;
					}
				}

				for ( int i = (int) len; i > 0; i-- ) {
					if ( endPos < 0 && !isblank( s[i - 1] ) ) {
						endPos = (int) i;
						break;
					}
				}

				if ( endPos > startPos ) {
					result = s.substr( startPos, endPos - startPos );
				}
			}

			return result;
		}

		/** Returns a copy of a string with all characters matching the cctype isblank predicate removed.*/
		static void TrimInPlace( string & s ) {
			size_t len = s.size();

			if ( len > 0 ) {
				int startPos = -1;
				int endPos = -1;

				for ( size_t i = 0; i < len; i++ ) {
					if ( startPos < 0 && !isblank( s[i] ) ) {
						startPos = (int) i;
						break;
					}
				}

				for ( int i = (int) len; i > 0; i-- ) {
					if ( endPos < 0 && !isblank( s[i - 1] ) ) {
						endPos = (int) i;
						break;
					}
				}

				if ( endPos > startPos ) {
					if ( (endPos < (int) len || startPos > 0) ) {
						s = s.substr( startPos, endPos - startPos );
					}
					// Else: string has no blank padding to remove.
				}
				else {
					// String contains nothing but spaces.
					s.clear();
				}
			}
		}

		template< typename T >
		static string Join( const T & collection, string delimiter = "," ) {
			string result;
			bool deja = false;

			for ( auto s : collection ) {
				if ( deja )
					result += delimiter;
				else
					deja = true;

				result += s;
			}

			return result;
		}


		static std::string Format( const char * format, ... ) {
			va_list args;
			va_start( args, format );
			char *buffer;
			std::size_t len;
			FILE *f = open_memstream( &buffer, &len );
			vfprintf( f, format, args );
			fclose( f );
			va_end( args );
			std::string result( buffer );
			free( buffer );
			return result;
		}
	};

}
