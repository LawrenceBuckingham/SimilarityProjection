#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <sstream>

#include "String.hpp"
#include "Array.hpp"
#include "Util.hpp"

using std::string;
using std::cerr;

using namespace QutBio;

// TODO: This clumsy attempt at a portable enum type can be used without too much change in JavaScript ports of this code. See if it can be made better.

namespace QutBio {

	/// <summary> Enumeration of the available 
	/// </summary>

	class EnumBase {
		friend ostream & operator<<( ostream &, const EnumBase & );

	private:
		std::string name;
		int value;

	protected:
		EnumBase( std::string name, int value ) : name( String::ToLowerCase( name ) ), value( value ) {}

	public:
		const string & ToString() {
			return name;
		}

		template<typename T>
		static T * Parse( const std::string & s, std::vector<EnumBase *> & values ) {
			string t = String::ToLowerCase( s );

			for ( size_t i = 0; i < values.size(); i++ ) {
				// cerr << "Checking " << values[i]->name << "\n";

				if ( t == values[i]->name ) {
					return (T*) values[i];
				}
			}

			ostringstream str;
			str << "Format Exception. Enumerated value '" << s << "' not recognised.";
			( cerr << str.str() << "\n" ).flush();
			throw new Exception( str.str(), __FILE__, __LINE__ );
		}

		bool operator == ( const EnumBase & other ) const {
			return this == &other ? true : ( this->value == other.value && this->name == other.name );
		}

		bool operator != ( const EnumBase & other ) const {
			return !( this == &other ? true : ( this->value == other.value && this->name == other.name ) );
		}

		int Value() { return value; }

		const string & Name() { return name; }
	};

	ostream & operator << ( ostream & out, const EnumBase & item ) {
		out << item.name;
		return out;
	}
}
