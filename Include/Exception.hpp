#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <stdexcept>

using string = std::string;
using runtime_error = std::runtime_error;

namespace QutBio {

	class Exception : public std::runtime_error {
	private:
		string file;
		int line;

	public:
		Exception( const string & message_, const string & file, int line )
			: runtime_error( message_ ), file( file ), line( line ) {}

		const string & File() { return file; };

		int Line() { return line; }
	};

	class FormatException : public Exception {
	public:
		FormatException( const string & message_, const string & file, int line )
			: Exception(message_, file, line ) {}
	};

	class KeyNotFoundException : public Exception {
	private:
		string key;

	public:
		KeyNotFoundException( const string & message_, const string & file, int line )
			: Exception( message_, file, line ), key( "Key not found" ) {}

		KeyNotFoundException( const string & message_, const string & key_, const string & file, int line )
			: Exception( message_, file, line ), key( key_ ) {}

		string & Key() { return key; }
	};

	class ArgumentException : public Exception {
	private:
		string arg;

	public:
		ArgumentException( const string & file, int line )
			: Exception( "Invalid argument", file, line ), arg( "Invalid argument" ) {}

		ArgumentException( const string & arg, const string & file, int line )
			: Exception( "Invalid argument: " + arg, file, line ), arg( arg ) {}

		string & Arg() { return arg; }
	};

	class NotImplementedException : public Exception {
	public:
		NotImplementedException( const string & file, int line )
			: Exception( "Not implemented.", file, line ) {}
	};
}

#define FileAndLine __FILE__, __LINE__

