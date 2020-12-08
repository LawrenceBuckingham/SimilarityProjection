#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Widget.H>
#include <FL/fl_draw.H>
#include <memory>
#include <string>
#include <vector>

using namespace std;

#include "LineSpec.hpp"

namespace LBGraph {

	/** Class hierarchy to render markers */
	class Marker
	{
	protected:
		/** Rendering pearameters for the marker outline. */
		LineSpec *line;

		/** The fill colour for the marker. */
		Fl_Color fillColour;

		/** The radius of the marker. Set this to negative value to suppress rendering of markers. */
		int radius;

	public:
		/**
			 * Constructor.
			 * @param line The line specification for drawing the outline.
			 * @param fillColour The colour to fill the shape (if that makes sense: some shapes are symbols.)
			 * @param radius The distance of the edge of the shape bounding box from the centre. If radius is zero, most marker implementations will probably draw a dot, and I expect they will generally only draw symmetrical shapes with an odd number of rows and columns in their bounding box.
			 */
		Marker(
			LineSpec *line = 0,
			Fl_Color fillColour = 0,
			int radius = -1
			//
		) : line(line ? line->Clone() : 0),
			fillColour(fillColour),
			radius(radius)
		{}

		/**
			 * Sets the line spec for the marker.
			 * @param line The new line specifier.
			 */
		void Line(LineSpec *line) { this->line = line; }

		/**
			 * Gets the current line specifier for the marker.
			 */
		LineSpec *Line() { return line; }

		/**
			 * Sets the fill colour for the marker.
			 * @param fillColour The new fill colour.
			 */
		void FillColour(Fl_Color fillColour) { this->fillColour = fillColour; }

		/**
			 * Gets the current line specifier for the marker.
			 */
		Fl_Color FillColour() { return fillColour; }

		/**
			 * Sets the fill colour for the marker.
			 * @param fillColour The new fill colour.
			 */
		void Radius(int radius) { this->radius = radius; }

		/**
			 * Gets the current line specifier for the marker.
			 */
		int Radius() { return radius; }

		/**
			 * Renders the marker at designated screen location.
			 * @param x The horizontal pixel location at which the marker will be centred.
			 * @param y The vertical pixel location at which the marker will be centered.
			 */
		virtual void Draw(int x, int y)
		{
			// default action: draw nothing.
		}

		/** Destructor.  */
		virtual ~Marker() {}

		/**
		 * Create a clone of the designated marker.
		 */
		virtual Marker *Clone() const
		{
			return new Marker(line, fillColour, radius);
		}
	};

	/** A Marker which is rendered as a 'x'. */
	class Cross : public Marker
	{
	public:
		/**
			 * Constructor.
			 * @param line The line specification for drawing the outline.
			 * @param radius The distance of the edge of the shape bounding box from the centre. If radius is zero, most marker implementations will probably draw a dot, and I expect they will generally only draw symmetrical shapes with an odd number of rows and columns in their bounding box.
			 */
		Cross(LineSpec *line, int radius = 3)
			: Marker(line, 0, radius)
		{}

		/**
			 * Renders a cross, centered on the screen location (x,y)
			 */
		void Draw(int x, int y) override
		{
			if (!line) return;

			int x0 = x - radius, x1 = x + radius + line->Thickness(), y0 = y - radius, y1 = y + radius + line->Thickness();
			fl_color(line->Colour());
			fl_line_style(line->Style(), line->Thickness());
			fl_line(x0, y0, x1, y1);
			fl_line(x0, y1, x1, y0);
			fl_line_style(0);
		}

		Marker * Clone() const override
		{
			return new Cross(line, radius);
		}
	};

	/** A Marker which is rendered as a '+'. */
	class Plus : public Marker
	{
	public:
		/**
			 * Constructor.
			 * @param line The line specification for drawing the outline.
			 * @param radius The distance of the edge of the shape bounding box from the centre. If radius is zero, most marker implementations will probably draw a dot, and I expect they will generally only draw symmetrical shapes with an odd number of rows and columns in their bounding box.
			 */
		Plus(LineSpec *line, int radius = 3)
			: Marker(line, 0, radius)
		{}

		/**
			 * Renders a '+', centered on the screen location (x,y)
			 */
		void Draw(int x, int y) override
		{
			if (!line) return;

			int x0 = x - radius, x1 = x + radius + line->Thickness(), y0 = y - radius, y1 = y + radius + line->Thickness();
			fl_color(line->Colour());
			fl_line_style(line->Style(), line->Thickness());
			fl_line(x, y0, x, y1);
			fl_line(x0, y, x1, y);
			fl_line_style(0);
		}

		Marker * Clone() const override
		{
			return new Plus(line, radius);
		}
	};

	class Circle : public Marker
	{
	public:
		/// Construct a Circle.
		///	@param line The address of a line specification which will be cloned into this object.
		///	@param fillColour The FLTK colour specifier for the filled body of the shape. If this is FL_INACTIVE_COLOR, then the obejct will not be filled.
		///	@param radius The number of pixels either side of the centre to be included in the shape. If radius is zero, a single pixel will be drawn (approximately, according to the detailed whim of the underlying FLTK rendering algorithm).
		Circle(LineSpec * line, Fl_Color fillColour, int radius = 3) : Marker(line, fillColour, radius) {}

		/// Draws the marker centered on the designated (x,y) screen location.
		/// @param x The absolute horizontal screen position at which the centre of the marker will appear.  
		/// @param y The absolute vertical screen position at which the centre of the marker will appear.  
		void Draw(int x, int y) override
		{
			if (fillColour != FL_INACTIVE_COLOR)
			{
				fl_color(fillColour);
				fl_begin_polygon();
				fl_circle(x, y, radius);
				fl_end_polygon();
			}

			if (line)
			{
				fl_color(line->Colour());
				fl_line_style(line->Style(), line->Thickness());
				fl_begin_line();
				fl_circle(x, y, radius);
				fl_end_line();
				fl_line_style(0);
			}
		}

		/// Creates a new dynamically allocated copy of this marker.
		///	@returns new Circle(line, fillColour, radius).
		Marker * Clone() const override
		{
			return new Circle(line, fillColour, radius);
		}
	};

	class Square : public Marker
	{
	public:
		/// Construct an axis-parallel Square marker.
		///	@param line The address of a line specification which will be cloned into this object.
		///	@param fillColour The FLTK colour specifier for the filled body of the shape. If this is FL_INACTIVE_COLOR, then the obejct will not be filled.
		///	@param radius The number of pixels either side of the centre to be included in the shape. If radius is zero, a single pixel will be drawn (approximately, according to the detailed whim of the underlying FLTK rendering algorithm).
		Square(LineSpec * line, Fl_Color fillColour, int radius = 3) : Marker(line, fillColour, radius) {}

		/// Draws the marker centered on the designated (x,y) screen location.
		/// @param x The absolute horizontal screen position at which the centre of the marker will appear.  
		/// @param y The absolute vertical screen position at which the centre of the marker will appear.  
		void Draw(int x, int y) override
		{
			if (fillColour != FL_INACTIVE_COLOR)
			{
				fl_color(fillColour);
				fl_rectf(x - radius, y - radius, radius * 2 + 1, radius * 2 + 1);
			}

			if (line)
			{
				fl_color(line->Colour());
				fl_line_style(line->Style(), line->Thickness());
				fl_loop(
					x - radius, y - radius,
					x + radius, y - radius,
					x + radius, y + radius,
					x - radius, y + radius
				);
				fl_line_style(0);
			}
		}

		/// Creates a new dynamically allocated copy of this marker.
		///	@returns new Circle(line, fillColour, radius).
		Marker * Clone() const override
		{
			return new Square(line, fillColour, radius);
		}
	};

	class Diamond : public Marker
	{};

	class Triangle : public Marker
	{};
}
