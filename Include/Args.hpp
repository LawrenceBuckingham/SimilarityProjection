#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "Exception.hpp"
#include "EnumBase.hpp"
#include "SimilarityMatrix.hpp"
#include "String.hpp"
#include "Util.hpp"

namespace QutBio {

	/// <summary> Helper class which parses command line arguments into
	///	a dictionary.
	/// </summary>
	class Args {
		friend std::ostream & operator<< ( std::ostream & str, Args args );

	private:
		unordered_map<string, vector<string>> arguments;

		string title;
		map<string, string> helpText;
		bool ok = true;

	public:

		Args( int argc, char ** argv ) : title( "" ) {
			helpText["help"] = "Generates this help text.";
			ParseArgs( argc, argv, arguments );
		}

		Args( size_t argc, const char ** argv ) : title( "" ) {
			helpText["help"] = "Generates this help text.";
			ParseArgs( (int) argc, (char**) argv, arguments );
		}

		Args( int argc, char ** argv, const string & title ) : title( title ) {
			helpText["help"] = "Generates this help text.";
			ParseArgs( argc, argv, arguments );
		}

		Args( size_t argc, const char ** argv, const string & title ) : title( title ) {
			helpText["help"] = "Generates this help text.";
			ParseArgs( (int) argc, (char**) argv, arguments );
		}

		string Title() {
			return this->title;
		}

		void SetTitle( const string & title ) {
			this->title = title;
		}

		void Help() {
			if ( !IsDefined( "help" ) ) return;

			for ( auto & s : helpText ) {
				if ( s.first.size() == 0 ) continue;

				cerr << "--" << s.first;
				auto parts = String::Split( s.second, "\n" );

				for ( auto &p : parts ) {
					cerr << "\n\t" << p;
				}

				cerr << "\n\n";
			}
		}

		bool Ok() {
			return ok;
		}

		void Fail() {
			ok = false;
		}

		void Reset() {
			ok = true;
		}

		bool Contains( string & key ) {
			string key_ = String::ToLowerCase( key );
			return arguments.find( key_ ) != arguments.end();
		}

		bool IsDefined( const char * key ) {
			string key_ = String::ToLowerCase( key );
			return arguments.find( key_ ) != arguments.end();
		}

		///// <summary> Safely gets the value of a single-valued argument. If the 
		/////		named argument is not present, returns undefined. Otherwise, returns 
		/////		the first element of the named list of arguments.
		///// </summary>
		///// <param name="arguments"></param>
		///// <param name="key"></param>
		///// <returns></returns>

		//bool Get(const string & key, string & result) {
		//	string key_(String::ToLowerCase(key));

		//	if (arguments.find(key_) != arguments.end()) {
		//		vector<string> & values(arguments[key_]);

		//		result = values.size() > 0 ? values[0] : "";
		//		return true;
		//	}

		//	return false;
		//}

		/// <summary> Safely gets the value of a single-valued argument. If the 
		///		named argument is not present, returns undefined. Otherwise, returns 
		///		the first element of the named list of arguments.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool GetImpl( const string & key, vector<string> & result, const string & help = "No help provided" ) {
			helpText[key] = help;

			string key_( String::ToLowerCase( key ) );

			if ( arguments.find( key_ ) != arguments.end() ) {
				vector<string> & values( arguments[key_] );
				result = values;
				return true;
			}

			return false;
		}

		/// <summary> Safely gets the value of a single-valued argument. If the 
		///		named argument is not present, returns false. Otherwise, returns 
		///		true and copies back by reference the first element of the 
		///		named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <param name="result"></param>
		/// <returns></returns>

		bool Get( const char * key, string & result, const string & help = "No help provided" ) {
			vector<string> vals;
			bool ok = GetImpl( string( key ), vals, help );

			if ( ok && vals.size() > 0 ) {
				result = vals[0];
				return true;
			}
			else {
				return false;
			}
		}

		/// <summary> Safely gets the value of a single-valued argument. If the 
		///		named argument is not present, returns false. Otherwise, returns 
		///		true and copies back by reference the first element of the 
		///		named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <param name="result"></param>
		/// <returns></returns>

		bool Get( const string & key, string & result, const string & help = "No help provided" ) {
			vector<string> vals;
			bool ok = GetImpl( key, vals, help );

			if ( ok && vals.size() > 0 ) {
				result = vals[0];
				return true;
			}
			else {
				return false;
			}
		}

		/// <summary> Replaces the value of an argument.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		void  Set( const string & key, const string & result ) {
			string key_( String::ToLowerCase( key ) );
			auto & results = arguments[key_];
			results.push_back( result );
		}

		/// <summary> Safely gets the value of a vector-valued argument. If the 
		///		named argument is not present, returns false. Otherwise, returns 
		///		true and passes back (by reference) the named list of values.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( string & key, vector<string> & result, const string & help = "No help provided" ) {
			vector<string> vals;
			result.clear();
			bool ok = GetImpl( key, vals, help );
			if ( ok ) result = vals;
			return ok;
		}

		/// <summary> Safely gets the value of a single-valued argument. If the 
		///		named argument is not present, returns undefined. Otherwise, returns 
		///		the first element of the named list of arguments.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key, vector<string> & result, const string & help = "No help provided" ) {
			string k( key );
			return Get( k, result, help );
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, double & result, const string & help = "No help provided" ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) && vals.size() > 0 ) {
				try {
					result = Double::Parse( vals[0] );
					return true;
				}
				catch ( Exception & ex ) {
					this->ok = false;
					cerr << "Invalid numeric data for argument " << key << ".\n";
					return false;
				}
			}
			else {
				return false;
			}
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="arguments"></param>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key, double & result, const string & help = "No help provided" ) {
			return Get( string( key ), result, help );
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as an int. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, int & result, const string & help = "No help provided" ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) && vals.size() > 0 ) {
				try {
					result = Int::Parse( vals[0] );
					return true;
				}
				catch ( Exception & ex ) {
					this->ok = false;
					cerr << "Invalid numeric data for argument " << key << ".\n";
					return false;
				}
			}
			else {
				return false;
			}
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as an int. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, uint & result, const string & help = "No help provided" ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) && vals.size() > 0 ) {
				try {
					result = Uint::Parse( vals[0] );
					return true;
				}
				catch ( Exception & ex ) {
					this->ok = false;
					cerr << "Invalid numeric data for argument " << key << ".\n";
					return false;
				}
			}
			else {
				return false;
			}
		}

#if DISTANCE_IS_NOT_INT
		/// <summary> Safely gets the value of a single-valued argument and parses it as an int. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, Distance & result, const string & help = "No help provided" ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) && vals.size() > 0 ) {
				try {
					result = Int::Parse( vals[0] );
					return true;
				}
				catch ( Exception &ex ) {
					this->ok = false;
					cerr << "Invalid numeric data for argument " << key << ".\n";
					return false;
				}
			}
			else {
				return false;
			}
		}
#endif

		/// <summary> Safely gets the value of a single-valued argument and parses it as a size_t. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, size_t & result, const string & help = "No help provided" ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) && vals.size() > 0 ) {
				try {
					result = Uint::Parse( vals[0] );
					return true;
				}
				catch ( Exception &ex ) {
					this->ok = false;
					cerr << "Invalid numeric data for argument " << key << ".\n";
					return false;
				}
			}
			else {
				return false;
			}
		}

		bool Get( const string & key, long & result, const string & help = "No help provided" ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) && vals.size() > 0 ) {
				try {
					istringstream str( vals[0] );
					str >> result;
					return true;
				}
				catch ( Exception &ex ) {
					this->ok = false;
					cerr << "Invalid numeric data for argument " << key << ".\n";
					return false;
				}
			}
			else {
				return false;
			}
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key, int & result, const string & help = "No help provided" ) {
			return Get( string( key ), result, help );
		}


		/// <summary> Safely gets the value of a single-valued argument and parses it as a boolean
		///		value. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key, bool & result, const string & help = "No help available." ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) && vals.size() > 0 ) {
				try {
					result = Bool::Parse( vals[0] );
					return true;
				}
				catch ( Exception &ex ) {
					this->ok = false;
					cerr << "Invalid numeric data for argument " << key << ".\n";
					return false;
				}
			}
			else {
				return false;
			}
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a boolean. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key, bool & result, const string & help = "No help available." ) {
			return Get( string( key ), result, help );
		}


		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns false. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const string & key ) {
			string keyLC( String::ToLowerCase( key ) );

			if ( arguments.find( keyLC ) != arguments.end() ) {
				auto values = arguments[keyLC];
				return values.size() > 0 ? String::ToLowerCase( values[0] ) == "true" : true;
			}
			else {
				return false;
			}
		}

		/// <summary> Safely gets the value of a single-valued argument and parses it as a double. 
		///		If the named argument is not present, returns undefined. Otherwise, returns the 
		///		first element of the named list of arguments.
		/// </summary>
		/// <param name="key"></param>
		/// <returns></returns>

		bool Get( const char * key ) {
			return Get( string( key ) );
		}

		void Show() {
			cout << (*this);
		}

		template<typename T>
		bool Get( const string & key, vector<EnumBase *> values, T *& result, const string & help = "No help available." ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) && vals.size() > 0 ) {
				try {
					result = EnumBase::Parse<T>( vals[0], values );
					return true;
				}
				catch ( Exception & ex ) {
					this->ok = false;
					cerr << "Invalid data for argument " << key << ".\n";
					return false;
				}
			}
			else {
				return false;
			}
		}

		template <typename T>
		bool Required( const string & key, vector<EnumBase *> values, T *& result, const string & help = "No help available." ) {
			if ( Get( key, values, result, help ) ) {
				return true;
			}
			else {
				ok = false;
				return false;
			}
		}

		template <typename T>
		bool Optional( const string & key, vector<EnumBase *> values, T *& result, const string & help = "No help available." ) {
			return Get( key, values, result, help );
		}

		template<typename T>
		bool Get( const char *key, unordered_set<T> & values, const string & help = "No help available." ) {
			vector<string> vals;

			if ( GetImpl( key, vals, help ) ) {
				for ( auto s : vals ) {
					istringstream str( s );
					T value;
					str >> value;
					values.insert( value );
				}

				return true;
			}
			else {
				return false;
			}
		}

		template <typename T>
		bool GetOptionalArgument( const char * name, T & arg ) {
			return IsDefined( name ) ? Get( name, arg ) : true;
		}

		template <typename T>
		bool GetOptionalArgument( const char * name, T & arg, Action showHelp ) {
			if ( IsDefined( name ) ) {
				if ( !Get( name, arg ) ) {
					cerr << "Unable to parse " << name << "." << endl;
					showHelp();
					return false;
				}
			}

			return true;
		}

		template <typename T>
		bool Get( const char * name, vector<T> & values, const string & help = "No help available." ) {
			vector<string> vals;

			if ( GetImpl( name, vals, help ) ) {
				for ( auto s : vals ) {
					istringstream str( s );
					T value;
					str >> value;
					values.push_back( value );
				}

				return true;
			}
			else {
				return false;
			}
		}

		/**
		 *	Gets a Similarity matrix based on an argument list, which will need to have
		 */

		bool Required( Alphabet *&alphabet, SimilarityMatrix *&matrix ) {
			helpText["matrixId"] = "Optional. The ID of a Blosum matrix.";
			helpText["matrixFile"] = "Optional. The name of a text file containing a matrix.";
			helpText["isCaseSensitive"] = "Optional, default = true. Is the similarity matrix case-sensitive.";

			int matrixId = 0;
			string id;

			if ( IsDefined( "matrixId" ) ) {
				if ( !Get( "matrixId", matrixId ) ) {
					cerr << ProgName() + ": error - argument 'matrixId' not valid.";
					ok = false;
					return false;
				}

				vector<int> matrices{ 35, 40, 45, 50, 62, 80, 100 };

				bool found = false;

				for ( auto x : matrices ) {
					if ( x == matrixId ) {
						found = true;
					}
				}

				if ( !found ) {
					cerr << ProgName() + ": error - matrix id not recognised.";
					ok = false;
					return false;
				}

				id = "matrixId " + Int::ToString( matrixId );
			}

			string matrixFile;
			DistanceType * distanceType = DistanceType::BlosumDistance();

			if ( IsDefined( "matrixFile" ) ) {
				Get( "matrixFile", matrixFile );
				distanceType = DistanceType::Custom();
				matrixId = -1;
				id = "--matrixFile '" + matrixFile + "'";
			}

			bool isCaseSensitive = false;

			if ( IsDefined( "isCaseSensitive" ) ) {
				if ( !Get( "isCaseSensitive", isCaseSensitive ) ) {
					cerr << ProgName() + ": Invalid data for argument 'isCaseSensitive'.";
					ok = false;
					return false;
				}

				id = id + " --isCaseSensitive " + (isCaseSensitive ? "'true'" : "'false'");
			}

			matrix = SimilarityMatrix::GetMatrix( alphabet, distanceType, matrixId, matrixFile );

			if ( !matrix ) {
				cerr << ProgName() + ": Unable to construct similarity matrix.\n";
				ok = false;
				return false;
			}

			matrix->id = id;
			alphabet = matrix->Alphabet();
			return true;
		}


		string ProgName() {
			string progName;
			Get( "", progName );
			return progName;
		}

		template <typename T>
		bool Required( T & value, const string & key, const string & help = "No help available." ) {
			helpText[key] = "Required. " + help;
			bool result = Get( key, value, help );

			if ( !result ) {
				ok = false;
				cerr << ProgName() << ": Required argument --" << key << " not found.\n";
			}

			return result;
		}

		template <typename T>
		bool Optional( T & value, const string & key, const string & help ) {
			helpText[key] = "Optional. " + help;
			bool result = Get( key, value, help );

			if ( !result ) {
				cerr << ProgName() << ": Optional argument --" << key << " not found. Using default value "
					<< value << ". \n";
			}

			return result;
		}

		/// <summary> Parses the supplied command line arguments and returns
		///		a unordered_map containing the results. For each named argument 
		///		(specified by a leading '-') a list of values is returned.
		///	<para>
		///		Arguments which appear before any named arguments are returned 
		///		in the entry with key = string.Empty.
		/// </para>
		/// </summary>
		/// <param name="args">
		///		A list of string containing the arguments.
		/// </param>
		/// <returns>
		///		A unordered_map with the arguments
		/// </returns>

	private:

		static void ParseArgs( int argc, char ** argv, unordered_map<string, vector<string>> & arguments ) {
			string currentKey;
			vector<string> * currentValues = &arguments[currentKey];

			for ( int i = 0; i < argc; i++ ) {
				string arg( argv[i] );

				if ( arg.size() >= 2 && arg[0] == '-' && arg[1] == '-' ) {
					currentKey = arg;
					currentKey.erase( currentKey.begin() );
					currentKey.erase( currentKey.begin() );
					currentValues = &arguments[String::ToLowerCase( currentKey )];
				}
				else {
					currentValues->push_back( arg );
				}
			}
		}
	};

	/// <summary> Echoes the command line arguments, parsed into lists.
	/// </summary>
	/// <param name="arguments"></param>

	ostream & operator<< ( ostream & str, Args args ) {
		bool deja( false );

		for ( auto arg : args.arguments ) {
			if ( deja ) str << " \\" << endl;

			deja = true;

			str << "--" << arg.first;

			for ( auto value : arg.second ) {
				str << " " << value;
			}
		}

		str << endl;

		return str;
	}

}
