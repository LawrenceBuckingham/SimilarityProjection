#pragma once
#include "IncludeAll.hpp"
#include <map>
#include <Alphabet.hpp>

using namespace std;

/**
 * Combo box to choose Alphabet type (DNA or RNA).
 */
class AlphabetChooser :
	public Fl_Choice {
	map<string, int> index;
	map<string, pAlphabet> stdAlphabets;
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
	AlphabetChooser( int left, int top, int width, int height ) :
		Fl_Choice( left, top, width, height, "DNA or RNA?" ) {
		stdAlphabets = Alphabets::StandardAlphabets();

		for ( auto & p : stdAlphabets ) {
			index[p.first] = add( p.first.c_str() );
		}

		Fl_Choice::value( 1 );
	}

	/**
		<summary>
			Returns the address of the currently selected Alphabet.
		</summary>
	 */
	Alphabet * value() const {
		// Want to override, but cannot because I'm trying to change the type
		auto selectedText = text();
		return selectedText == nullptr ? nullptr : stdAlphabets.at(selectedText);
	}

	/**
	 *	Selects a designated Alphabet.
	 *	@param val The address of the Alphabet that will be pre-selected.
	 *	@pre val in ran(Alphabets::StandardAlphabets())
	 *	@post value() == val
	 */
	void value( Alphabet * val ) {
		for ( auto p : stdAlphabets ) {
			if ( p.second == val ) {
				Fl_Choice::value( index[p.first] );
			}
		}
	}

	/**
	 *	Sets the value of this control based on the alphabet name.
	 *	@param	name The name of a standard alphabet.
	 *	@pre	name in dom(Alphabets::StandardAlphabets())
	 *	@post	value() == Alphabets::StandardAlphabets()[name]
	 */
	void value( const string & name ) {
		Fl_Choice::value( index[name] );
	}
};
