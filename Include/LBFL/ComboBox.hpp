#pragma once

#include <cstdarg>
#include <map>
#include <string>
#include <vector>
#include <FL/Fl_Choice.H>

using namespace std;

/**
 * Combo box to choose Alphabet type (DNA or RNA).
 */
class ComboBox : public Fl_Choice {
	vector<string> values;
	map<string, int> index;

public:
	/**
		<summary>
			DnaOrRna constructor.
		</summary>
		<param name="left">
			The nominal position of the left edge of the control.
		</param>
		<param name="top">
			The nominal position of the left edge of the control.
		</param>
		<param name="width">
			The desired width of the control.
		</param>
		<param name="height">
			The desired width of the control.
		</param>
	 */
	ComboBox( int left, int top, int width, int height ) :
		Fl_Choice( left, top, width, height ) {}

	/**
	 *
	 */
	ComboBox & SetChoices( const char * s, ... ) {
		Fl_Choice::clear();
		values.clear();
		index.clear();
		values.emplace_back( s );
		Fl_Choice::add(values.back().c_str());
		index[s] = 0;

		va_list args;
		va_start(args, s);

		for ( auto val = va_arg( args, const char * ); val; val = va_arg( args, const char * ) ) {
			index[val] = values.size();
			values.emplace_back( val );
			Fl_Choice::add(values.back().c_str());
		}

		va_end(args);

		Fl_Choice::value( 0 );

		return *this;
	}
	
	/// Destructor
	virtual ~ComboBox() {}

	/**
	 *	Sets the value of this control based on the alphabet name.
	 *	@param	name The name of a standard alphabet.
	 *	@pre	name in dom(values)
	 *	@post	value() == values[name]
	 */
	void value( const string & name ) {
		Fl_Choice::value( index[name] );
	}

	int value () const {
		return Fl_Choice::value();
	}

	void value ( int val ) {
		Fl_Choice::value(val);
	}
};
