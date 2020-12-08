#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Widget.H>
#include <FL/fl_draw.H>
#include <memory>
#include <string>
#include <vector>

using namespace std;

#include "FontSpec.hpp"
#include "Label.hpp"
#include "LineSpec.hpp"
#include "Marker.hpp"
#include "Projection.hpp"

namespace LBGraph
{
class Point
{
	Label *label;	//< Label definition for this point.
	Marker *marker;  //< Point custom marker. Leave this as null if you want to use a default marker supplied by Series::Draw().
public:
	double x;  //< Position of point in the x-axis frame of reference.
	double y;  //< Position of point in the y-axis frame of reference.
	double z;  //< Third channel.

	/**
	 * Construct a new Point.
	 */
	Point(
		double x,	 //< x value.
		double y,	 //< y value.
		double z = 0  //< optional z value (a third channel if desired).
		) : x(x),
			y(y),
			z(z),
			label(0),
			marker(0)
	{
	}

	Point(const Point &other) : x(other.x),
								y(other.y),
								z(other.z),
								label(other.label ? other.label->Clone() : 0),
								marker(other.marker ? other.marker->Clone() : 0)
	{
	}

	Point &operator=(const Point &other)
	{
		x = other.x;
		y = other.y;
		z = other.z;
		label = other.label ? other.label->Clone() : 0;
		marker = other.marker ? other.marker->Clone() : 0;
		return *this;
	}

	~Point()
	{
		delete label;
		delete marker;
	}

	Point & SetLabel(
		const string &text,	//< the text that will be displayed.
		int offsetX = 0,	   //< Horizontal pixel displacement from screen-location of point and label anchor point.
		int offsetY = 0,	   //< Vertical pixel displacement from screen-location of point and label anchor point.
		double anchorX = 0.5,  //< Relative horizontal location of the label anchor point on the label bounding box. 0 == left edge of label bounding box, 1 == right edge of label bounding box.
		double anchorY = 0.5,  //< Relative vertical location of the label anchor point on the label bounding box. 0 == top edge of label bounding box, 1 == bottom edge of label bounding box.
		FontSpec *font = 0	 //< Point custom font.
	)
	{
		label = new LBGraph::Label(text, offsetX, offsetY, anchorX, anchorY, font);
		return *this;
	}

	Point & SetMarker(const Marker *prototype)
	{
		marker = prototype->Clone();
		return *this;
	}

	/**
		 * Renders the marker and/or label of the point at the designated screen location.
		 * @param x_ horizontal pixel location at which the centre of the marker would appear.
		 * @param y_ vertical pixel location at which the centre of the marker would appear.
		 */
	void Draw(
		Projection projection, 
		LBGraph::Marker *defaultMarker = 0, 
		FontSpec *defaultFont = 0
	) {
		int x, y;

		if (projection(this->x, this->y, x, y))
		{
			if (marker)
			{
				marker->Draw(x, y);
			}
			else if (defaultMarker)
			{
				defaultMarker->Draw(x, y);
			}

			if (label)
			{
				label->Draw(x, y, defaultFont);
			}
		}
	}
};
}  // namespace HBGraph
