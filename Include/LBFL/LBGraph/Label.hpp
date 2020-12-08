#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Widget.H>
#include <FL/fl_draw.H>
#include <memory>
#include <string>
#include <vector>

using namespace std;

#include "FontSpec.hpp"
#include "LineSpec.hpp"
#include "Marker.hpp"

namespace LBGraph
{
class Label
{
	string text;	 //< Text to display near point.
	int offsetX;	 //< Horizontal pixel displacement from screen-location of point and label anchor point.
	int offsetY;	 //< Vertical pixel displacement from screen-location of point and label anchor point.
	double anchorX;  //< Relative horizontal location of the label anchor point on the label bounding box. 0 == left edge of label bounding box, 1 == right edge of label bounding box.
	double anchorY;  //< Relative vertical location of the label anchor point on the label bounding box. 0 == top edge of label bounding box, 1 == bottom edge of label bounding box.
	FontSpec *font;  //< Private, owned, pointer to font specifier, or NULL is a default is to be used.

public:
	/**
	 *	Construct a new label. 
	 */
	Label(
		const string &text,	//< the text that will be displayed.
		int offsetX = 0,	   //< Horizontal pixel displacement from screen-location of point and label anchor point.
		int offsetY = 0,	   //< Vertical pixel displacement from screen-location of point and label anchor point.
		double anchorX = 0.5,  //< Relative horizontal location of the label anchor point on the label bounding box. 0 == left edge of label bounding box, 1 == right edge of label bounding box.
		double anchorY = 0.5,  //< Relative vertical location of the label anchor point on the label bounding box. 0 == top edge of label bounding box, 1 == bottom edge of label bounding box.
		FontSpec *font = 0	 //< Point custom font. This will be cloned if non-zero.

		) : text(text),
			offsetX(offsetX),
			offsetY(offsetY),
			anchorX(anchorX),
			anchorY(anchorY)
	{
		this->font = font ? font->Clone() : 0;
	}

	Label(const Label &other)
	{
		text = other.text;
		offsetX = other.offsetX;
		offsetY = other.offsetY;
		anchorX = other.anchorX;
		anchorY = other.anchorY;
		font = other.font ? other.font->Clone() : 0;
	}

	Label & operator=(const Label &other)
	{
		text = other.text;
		offsetX = other.offsetX;
		offsetY = other.offsetY;
		anchorX = other.anchorX;
		anchorY = other.anchorY;
		font = other.font ? other.font->Clone() : 0;
		return *this; 
	}

	virtual ~Label()
	{
		delete font;
	}

	void Draw(int x, int y, FontSpec *defaultFont = 0)
	{
		FontSpec *font = this->font;

		if (!font) font = defaultFont;

		if (font) fl_font(font->Family(), font->Size());

		int w = 0, h = 0;
		fl_measure(text.c_str(), w, h);

		int labelX = x + offsetX - anchorX * w;
		int labelY = y + offsetY - anchorY * h - fl_descent() + fl_height();
		fl_draw(text.c_str(), labelX, labelY);
	}

	Label *Clone() {
		return new Label(*this);
	}
};
}  // namespace HBGraph