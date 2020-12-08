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

/** Encapsulate a font specifier. */
class FontSpec
{
protected:
	Fl_Font family = 0;  //< The font family. These are enumerated in <FL/Enumerations.H>. This default is a Sans Serif system font.
	int size = 12;		 //< The pixel (not point) height of the font.
public:
	FontSpec( 
		Fl_Font family = 0, //< The font family. These are enumerated in <FL/Enumerations.H>. This default is a Sans Serif system font.
		int size = 12 		//< The pixel (not point) height of the font.
	) : family(family), size(size)
	{}

	void Family(Fl_Font family) { this->family = family; }

	void Size(int size) { this->size = size; }

	Fl_Font Family() { return family; }

	int Size() { return size; }

	/** Destructor. */
	virtual ~FontSpec() {}

	/** Gets a default font specifier. */
	static FontSpec *Default()
	{
		static FontSpec *instance(new FontSpec());
		return instance;
	}

	virtual FontSpec *Clone() const {
		return new FontSpec(this->family, this->size);
	};
};
}