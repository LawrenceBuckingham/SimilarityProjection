#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Widget.H>
#include <FL/fl_draw.H>
#include <memory>
#include <string>
#include <vector>

using namespace std;

namespace LBGraph
{
class LineSpec
{
protected:
	/** The border color for the point markers. */
	Fl_Color colour = FL_BLACK;

	/**
		 *	The thickness of the marker border line, measured in pixels. If zero, a reasonable
		 *	system-specific default thickness is used.
		 */
	int thickness = 0;

	/** The FLTK line style of the marker border line. Copied from FLTK documentation: A bitmask which is a bitwise-OR of a line style, a cap style, and a join style. If you don't specify a dash type you will get a solid line. If you don't specify a cap or join type you will get a system-defined default of whatever value is fastest. */
	int style = FL_SOLID;

	/** Byte pattern for dashes, encoded in a ASCIIZ string. */
	string dashes = "";

public:
	LineSpec(
		Fl_Color colour = FL_BLACK,
		int thickness = 0,
		int style = FL_SOLID,
		const string &dashes = "") : colour(colour), thickness(thickness), style(style), dashes(dashes)
	{
	}

	/** Destructor. */
	virtual ~LineSpec() {}

	void Colour(Fl_Color colour) { this->colour = colour; }
	Fl_Color Colour() { return colour; }

	void Thickness(int thickness) { this->thickness = thickness; }
	int Thickness() { return this->thickness; }

	void Style(int style) { this->style = style; }
	int Style() { return this->style; }

	string &Dashes() { return this->dashes; }

	/**
		 * Gets a useable default line specifier. Colour is black, thickness is zero (yielding a FLTK default thickness), style is solid and dash pattern is an empty sequence.
		 * @returns a shared pointer to a global default line specifier.
		 */
	static LineSpec *Default()
	{
		static LineSpec *instance(new LineSpec());
		return instance;
	}

	/**
	 * Gets a copy of this object.
	 * @returns A newly created clone of the present object.
	 */
	LineSpec *Clone() const {
		return new LineSpec(this->colour, this->thickness, this->style, this->dashes);
	}
};

}  // namespace HBGraph
