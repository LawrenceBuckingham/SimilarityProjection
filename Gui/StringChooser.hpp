#pragma once
#include "IncludeAll.hpp"
#include <map>
#include <Alphabet.hpp>

using namespace std;

/**
 * Combo box to choose Alphabet type (DNA or RNA).
 */
class StringChooser :
	public Fl_Choice {
	map<string, int> index;
	vector<string> values;
public:
	StringChooser( int left, int top, int width, int height, const char * label = 0 ) :
		Fl_Choice( left, top, width, height, label )
		{}

	void SetValues ( const vector<string> & values ) {
		this->values = values;
		
		for ( auto & p : this->values ) {
			index[p] = add( p.c_str() );
		}

		Fl_Choice::value( 0 );
	}

	string value() const {
		// Want to override, but cannot because I'm trying to change the type
		auto selectedText = text();
		return selectedText == nullptr ? "" : selectedText;
	}

	void value( const string & name ) {
		Fl_Choice::value( index[name] );
	}
};
