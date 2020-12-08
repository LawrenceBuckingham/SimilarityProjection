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
#include "Point.hpp"
#include "Projection.hpp"

namespace LBGraph {
	class Series {
	protected:
		vector<Point *> data;
		Marker *marker = 0;  //< Default marker for points in series.
		FontSpec *font = 0;  //< Default typeface for labels in series.
		LineSpec *line = 0;  //< Line style for scatter plots.
		string name;	 //< Series name.
	public:
		Series(
			Marker *marker = 0,
			LineSpec *line = 0,
			FontSpec *font = 0,
			const string &name = ""
			//
		) : name( name ) //
		{
			SetMarker( marker );
			SetLine( line );
			SetFont( font );
		}

		virtual ~Series() {
			Clear();
			delete line;
			delete font;
			delete marker;
		}

		Point *Add( double x, double y, double z = 0 ) {
			Point *p = new Point( x, y, z );
			data.push_back( p );
			return p;
		}

		void Clear() {
			for ( auto p : data ) {
				delete p;
			}
			data.clear();
		}

		/** Gets a reference to the points in the series. */
		const vector<Point *> & Data() const {
			return data;
		}

		/**
		 *  Overwrites the current marker symbol of the series with a clone of the supplied argument.
		 * @param marker The address of a marker that (if not null) will be cloned into this series. If this is null, the existing marker (if any) will be deleted and overwritten with null.
		 */
		void SetMarker( const Marker *marker ) {
			if ( marker == this->marker ) return;

			delete this->marker;
			this->marker = marker ? marker->Clone() : 0;
		}

		/**
		 *  Overwrites the current font of the series with a clone of the supplied argument.
		 * @param font The address of a font specifier that (if not null) will be cloned into this series. If this is null, the existing font specifier (if any) will be deleted and overwritten with null.
		 */
		void SetFont( const FontSpec *font ) {
			if ( font == this->font ) return;

			delete this->font;
			this->font = font ? font->Clone() : 0;
		}

		/**
		 *  Overwrites the current line style of the series with a clone of the supplied argument.
		 * @param font The address of a line specifier that (if not null) will be cloned into this series. If this is null, the existing line specifier (if any) will be deleted and overwritten with null.
		 */

		void SetLine( const LineSpec *line ) {
			if ( line == this->line ) return;

			delete this->line;
			this->line = line ? line->Clone() : 0;
		}

		/**
		 *	Gets the address of the marker (if any). Returns null if there is no registered marker.
		 *	@returns The address of the marker (if any). Returns null if there is no registered marker.
		 */
		LBGraph::Marker *GetMarker() { return marker; }

		/**
		 *	Gets the address of the marker (if any). Returns null if there is no registered marker.
		 *	@returns The address of the marker (if any). Returns null if there is no registered marker.
		 */
		FontSpec *Font() { return font; }

		/**
		 *	Gets the address of the marker (if any). Returns null if there is no registered marker.
		 *	@returns The address of the marker (if any). Returns null if there is no registered marker.
		 */
		LineSpec *Line() { return line; }

		/**
		 * Gets a read-write reference to the name of the series.
		 * @returns A reference to the name field of this object.
		 */

		string &Name() { return name; }

		/**
			 * Draws the points and lines making up the series.
			 */
		virtual void Draw( Projection project ) {
			if ( line && line->Style() != -1 ) {
				for ( uint i = 1; i < data.size(); i++ ) {
					DrawLine( data[i - 1]->x, data[i - 1]->y, data[i]->x, data[i]->y,
						project, line->Colour(), line->Style(), line->Thickness(),
						line->Dashes().c_str() );
				}
			}

			for ( uint i = 0; i < data.size(); i++ ) {
				data[i]->Draw( project, marker, font );
			}
		}

		void DrawLine(
			double x0,
			double y0,
			double x1,
			double y1,
			Projection project,
			Fl_Color colour = 0,
			int lineStyle = 0,
			int thickness = 0,
			const string &dashes = "" ) {
			int i0, j0, i1, j1;

			if ( project( x0, y0, i0, j0 ) && project( x1, y1, i1, j1 ) ) {
				auto oldColour = fl_color();
				fl_color( colour );
				fl_line_style( lineStyle, thickness, (char *) dashes.c_str() );
				fl_line( i0, j0, i1, j1 );
				fl_color( oldColour );
				fl_line_style( 0 );
			}
		}

		pair<int, double> NearestTo( double x, double y ) {
			pair<int, double> res{ -1, numeric_limits<double>::max() };

			for ( int i = 0; i < (int) data.size(); i++ ) {
				auto p = data[i];
				double d2 = (x - p->x)*(x - p->x) + (y - p->y)*(y - p->y);

				if ( d2 < res.second ) {
					res.first = i;
					res.second = d2;
				}
			}

			res.second = sqrt( res.second );

			return res;
		}
	};
}  // namespace HBGraph
