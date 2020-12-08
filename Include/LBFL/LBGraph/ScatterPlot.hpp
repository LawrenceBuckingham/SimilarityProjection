#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Widget.H>
#include <FL/fl_draw.H>
#include <memory>
#include <string>
#include <vector>

using namespace std;

#include "Axis.hpp"
#include "FontSpec.hpp"
#include "Label.hpp"
#include "LineSpec.hpp"
#include "Marker.hpp"
#include "Point.hpp"
#include "Series.hpp"
#include "IEventHandler.hpp"
#include "Rectangle.hpp"
#include <cmath>

namespace LBGraph {
	/** A Graph is a widget that can plot series of points.*/
	class ScatterPlot :
		public Fl_Widget,
		public virtual IEventHandler {

	public:
		struct MouseEvent {
			ScatterPlot * source;
			double x;
			double y;
			int eventCode;
			int keyCode;
		};

		struct MouseHandler {
			virtual ~MouseHandler() {}
			virtual bool Handle( const MouseEvent & event ) = 0;
		};

	private:
		LinearAxis xAxis_, yAxis_;
		Axis *xAxis;
		Axis *yAxis;
		double xCrossesY = 0;
		double yCrossesX = 0;
		vector<Series *> series;
		int marginLeft, marginTop, marginRight, marginBottom;
		Fl_Color bgColour, fgColour;
		Series hTics;
		Series vTics;
		vector<IEventHandler *> eventHandlers;
		vector<MouseHandler *> mouseHandlers;

	public:
		ScatterPlot( int x = 0, int y = 0, int w = 0, int h = 0, const char *label = 0 )
			: Fl_Widget( x, y, w, h, label ),
			bgColour( color() ),
			fgColour( FL_BLACK ),
			xAxis_( 0, 1 ),
			yAxis_( 1, 0 ),
			xAxis( &xAxis_ ),
			yAxis( &yAxis_ )
			//
		{
			AddEventHandler( this );
		}

		virtual ~ScatterPlot() {}

		int handle( int eventCode ) override {
			int sx = Fl::event_x();
			int sy = Fl::event_y();

			int xx = x(), yy = y(), ww = w(), hh = h();

			if ( sx < xx || sy < yy || sx >= xx + ww | sy >= yy + hh ) return 0;

			for ( auto h : eventHandlers ) {
				if ( h->HandleEvent( this, eventCode ) ) return 1;
			}

			return 0;
		}

		ScatterPlot & AddEventHandler( IEventHandler * handler ) {
			eventHandlers.push_back( handler );
			return *this;
		}

		/**
		 * Gets the address of the current x-axis object. This initially refers
		 *	to a default LinearAxis, but that can be replaced via SetXAxis.
		 * @returns The address of the current x-axis object.
		 */
		Axis *XAxis() { return xAxis; }

		/**
		 * Set the new x-axis, or reset it back to the internal default.
		 * @param axis The address of an Axis object which will project x-coordinates.
		 *		A non-zero pointer is saved and used as-is, so it needs to refer
		 *		to an object having life-span at least as great as that of the
		 *		chart. If the value is 0, then the internal default (LinearAxis)
		 *		object is restored.
		 */
		ScatterPlot &  SetXAxis( Axis *axis = 0 ) {
			xAxis = axis ? axis : &xAxis_;
			return *this;
		}

		/**
		 * Gets the address of the current y-axis object. This initially refers
		 *	to a default LinearAxis, but that can be replaced via SetYAxis.
		 * @returns The address of the current y-axis object.
		 */
		Axis *YAxis() { return yAxis; }

		/**
		 * Set the new yAxis, remembering that y-axes need to have the limits reversed!!!
		 * @param axis The address of an Axis object which will project y-coordinates.
		 *		A non-zero pointer is saved and used as-is, so it needs to refer
		 *		to an object having life-span at least as great as that of the
		 *		chart. If the value is 0, then the internal default (LinearAxis)
		 *		object is restored.
		 */
		ScatterPlot &  SetYAxis( Axis *axis = 0 ) {
			yAxis = axis ? axis : &yAxis_;
			return *this;
		}

		/**
		 * Sets the current colour for axis and label drawing.
		 * @param foregroundColour new colour for axis and label drawing.
		 */
		ScatterPlot &  Colour( Fl_Color foregroundColour ) {
			fgColour = foregroundColour;
			return *this;
		}

		/**
		 * Gets the current colour for axis and label drawing.
		 * @returns the current foreground colour.
		 */
		Fl_Color Colour() {
			return fgColour;
		}

		/**
		 * Sets the colour for drawing the chart background.
		 * @param fillColour The new colour to use to fill the chart area.
		 */
		ScatterPlot & SetFillColour( Fl_Color fillColour ) {
			bgColour = fillColour;
			return *this;
		}

		/**
		 * Gets the current background (fill) colour.
		 * @returns the current background colour of the chart.
		 */
		Fl_Color FillColour() {
			return bgColour;
		}

		/**
		 * Adds the designated series to the collection. Simply adds the supplied address to the list of series, so any changes in the original series will be reflected in the chart when it is drawn.
			Equivalent to <code>Series().push_back(newSeries)</code>
		 * @param newSeries The address of a new series to add.
		 */
		ScatterPlot & Add( Series *newSeries ) {
			series.push_back( newSeries );
			return *this;
		}

		/**
		 *	Adds a mouse event handler to the list of objects that will be notified when mouse events are detected on this drawing surface.
		 *	@param handler The address of a MouseHandler.
		 *	@pre handler != 0.
		 *	@post handler in mouseHandlers.
		 */
		ScatterPlot & AddMouseHandler( MouseHandler * handler ) {
			mouseHandlers.push_back( handler );
			return *this;
		}

		/// Gets a reference to the list of stored data series.
		///	@returns this->series.
		vector<Series *> & DataSeries() {
			return this->series;
		}

		/**
		 * Gets a reference to the series used to render tick marks on the x-axis.
		 * @returns A reference to the hTics series.
		 */
		Series *HTics() {
			return &hTics;
		}

		/**
		 * Gets a reference to the series used to render tick marks on the x-axis.
		 * @returns A reference to the hTics series.
		 */
		Series *VTics() {
			return &vTics;
		}

		/**
			 * Gets the locations on each axis where they cross the other.
			 *
			 * @param yCrossesX reference variable that will be updated with the x-value at which the
			 * 	 		y-axis crosses the x-axis.
			 * @param xCrossesY reference variable that will be updated with the y-value at which the
			 * 	 		x-axis crosses the y-axis.
			 */
		ScatterPlot & AxesCross( double &yCrossesX, double &xCrossesY ) {
			yCrossesX = this->yCrossesX;
			xCrossesY = this->xCrossesY;
			return *this;
		}

		/**
			 * Set the locations on each axis where they cross the other.
			 *
			 * @param yCrossesX reference variable that will be updated with the x-value at which the
			 * 	 		y-axis crosses the x-axis.
			 * @param xCrossesY reference variable that will be updated with the y-value at which the
			 * 	 		x-axis crosses the y-axis.
			 */
		ScatterPlot &  SetAxesCross( double yCrossesX, double xCrossesY ) {
			this->yCrossesX = yCrossesX;
			this->xCrossesY = xCrossesY;
			return *this;
		}

		ScatterPlot & SetAxisBounds( double xLeft, double xRight, double yTop, double yBottom ) {
			xAxis->SetBounds( xLeft, xRight );
			yAxis->SetBounds( yTop, yBottom );
			return *this;
		}

		ScatterPlot &  SetMargin( int left, int top, int right, int bottom ) {
			this->marginLeft = left;
			this->marginTop = top;
			this->marginRight = right;
			this->marginBottom = bottom;
			return *this;
		}

		/**
			 * Removes all data and dsiplays an empty chart.
			 */
		ScatterPlot &  Clear() {
			series.clear();
			return *this;
		}

		Rectangle<int> PlotArea() {
			int left = x(), top = y(), right = left + w() - 1, bottom = top + h() - 1;
			Rectangle<int> result{ left + marginLeft, top + marginTop, right - marginRight, bottom - marginBottom };
			return result;
		}

		/**
			 * Helper function to draw a line between the designated endpoints using the current axes.
			 * This will not be remembered later unles the line is rendered as part of a series.
			 * @param x0 the position of the first endpoint on the x-axis.
			 * @param y0 the position of the first endpoint on the y-axis.
			 * @param x1 the position of the second endpoint on the x-axis.
			 * @param y1 the position of the second endpoint on the y-axis.
			 * @param colour the colour to render the line
			 * @param linestyle the style for the line. This needs to be one of the FLTK line styles.
			 * @param thickness the thickness to the line.
			 */
		void DrawLine(
			double x0,
			double y0,
			double x1,
			double y1,
			Fl_Color colour = 0,
			int lineStyle = 0,
			int thickness = 0,
			const string &dashes = "" ) {
			auto oldColour = fl_color();
			int left = x(), top = y(), right = left + w() - 1, bottom = top + h() - 1;
			int i0 = round( xAxis->ToScreen( x0, left + marginLeft, right - marginRight ) );
			int i1 = round( xAxis->ToScreen( x1, left + marginLeft, right - marginRight ) );
			int j0 = round( yAxis->ToScreen( y0, top + marginTop, bottom - marginBottom ) );
			int j1 = round( yAxis->ToScreen( y1, top + marginTop, bottom - marginBottom ) );
			fl_color( colour );
			fl_line_style( lineStyle, thickness, (char *) dashes.c_str() );
			fl_line( i0, j0, i1, j1 );
			fl_color( oldColour );
			fl_line_style( 0 );
		}

		function<bool( double wx, double wy, int &sx, int &sy )> Projection() {
			return [this]( double wx, double wy, int &sx, int &sy ) {
				int left = x(), top = y(), right = left + w() - 1, bottom = top + h() - 1;

				if ( xAxis->Min() <= wx && wx <= xAxis->Max() && yAxis->Max() <= wy && wy <= yAxis->Min() ) {
					sx = xAxis->ToScreen( wx, left + marginLeft, right - marginRight );
					sy = yAxis->ToScreen( wy, top + marginTop, bottom - marginBottom );
					return true;
				}
				else {
					return false;
				}
			};
		}

		/**
			 * Renders the chart.
			 */
		void draw() override {
			draw_box();
			Fl_Boxtype b = box();
			int xx = x() + Fl::box_dx( b );  // was 9 instead of dx...
			int yy = y() + Fl::box_dy( b );
			int ww = w() - Fl::box_dw( b );
			int hh = h() - Fl::box_dh( b );
			fl_push_clip( xx, yy, ww, hh );

			fl_color( bgColour );
			fl_rectf( xx, yy, ww, hh );

			if ( xAxis && yAxis ) {
				//fprintf( stderr, "xAxis: (%f,%f)\n", xAxis->Min(), xAxis->Max() );
				//fprintf( stderr, "yAxis: (%f,%f)\n", yAxis->Min(), yAxis->Max() );

				for ( auto ser : series ) {
					ser->Draw( Projection() );
				}

				if ( xAxis->IsVisible() ) {
					DrawLine( xAxis->Min(), xCrossesY, xAxis->Max(), xCrossesY );
				}

				if ( yAxis->IsVisible() ) {
					DrawLine( yCrossesX, yAxis->Min(), yCrossesX, yAxis->Max() );
				}

				hTics.Draw( Projection() );
				vTics.Draw( Projection() );
			}

			draw_label();
			fl_pop_clip();
		}

		struct PlottedPoint {
			Series * series;
			int index;
			double distance;
		};

		PlottedPoint NearestTo( double x, double y ) {
			PlottedPoint res{ 0, -1, numeric_limits<double>::max() };

			for ( auto ser : series ) {
				auto r = ser->NearestTo( x, y );

				if ( r.second < res.distance ) {
					res.series = ser;
					res.index = r.first;
					res.distance = r.second;
				}
			}

			return res;
		}

		bool HandleEvent( Fl_Widget *src, int eventCode ) override {
			// fprintf(stderr, "Received event %d\n", eventCode);

			if ( mouseHandlers.size() > 0 ) {
				int sx = Fl::event_x();
				int sy = Fl::event_y();

				int xx = x(), yy = y(), ww = w(), hh = h();

				if ( sx < xx || sy < yy || sx >= xx + ww | sy >= yy + hh ) return false;

				auto plotArea = PlotArea();
				double x = xAxis->ToWorld( sx, plotArea.left, plotArea.right );
				double y = yAxis->ToWorld( sy, plotArea.top, plotArea.bottom );

				MouseEvent me{ this, x, y, eventCode, Fl::event_key() };

				bool handled = false;

				for ( auto h : mouseHandlers ) {
					handled = handled || h->Handle( me );
				}

				return handled;
			}
			else {
				return false;
			}
		}


		ScatterPlot & SetTickMarks( int numHTicks, const char * hFmt, int numVTicks, const char * vFmt ) {
			LineSpec tickLine;
			Plus tickMark( &tickLine, 5 );
			FontSpec mono12{ FL_SCREEN, 12 };

			if ( hTics.GetMarker() == 0 ) {
				hTics.SetMarker( &tickMark );
			}

			hTics.Clear();

			bool isLogX = (0 != dynamic_cast<LogarithmicAxis *>(xAxis));

			double xMax = std::max( xAxis->Max(), xAxis->Min() );
			double xMin = std::min( xAxis->Max(), xAxis->Min() );

			for ( int i = 0; i <= numHTicks; i++ ) {
				double x = isLogX ? xMin * pow( 10, i ) : xMin + i * (xMax - xMin) / numHTicks;
				Point *p = hTics.Add( x, 0 );
				p->SetLabel( String::Format( hFmt, x ), 0, 7, 0.5, 0.0, &mono12 );
			}

			if ( vTics.GetMarker() == 0 ) {
				vTics.SetMarker( &tickMark );
			}

			vTics.Clear();

			// TODO: Move this into the Axis subclasses.
			bool isLogY = (0 != dynamic_cast<LogarithmicAxis *>(yAxis));

			double yMax = std::max( yAxis->Max(), yAxis->Min() );
			double yMin = std::min( yAxis->Max(), yAxis->Min() );

			for ( int i = 0; i <= numVTicks; i++ ) {
				double y = isLogY ? yMin * pow( 10, i ) : i * (yMax - yMin) / numVTicks;
				Point *p = vTics.Add( 0, y );
				p->SetLabel( String::Format( vFmt, y ), -7, 0, 1.0, 0.5, &mono12 );
			}

			return *this;
		}

		void SaveToFile( const string & fileName ) {
			uint maxRows = 0;
			auto & series = this->series;
			auto S = series.size();

			for ( const auto & s : series ) {
				uint n = s->Data().size();
				if ( n > maxRows ) maxRows = n;
			}

			if ( maxRows == 0 ) return;

			ofstream f( fileName );
			CsvWriter w( f );

			for ( uint i = 0; i < S; i++ ) {
				auto s = series[i];
				auto name = s->Name().c_str();
				string xName = String::Format( "%s_x", name );
				string yName = String::Format( "%s_y", name );
				w << xName << yName;
			}

			w.Ln();

			for ( uint row = 0; row < maxRows; row++ ) {
				for ( uint col = 0; col < S; col++ ) {
					auto s = series[col];
					auto & data = s->Data();

					if ( row < data.size() ) {
						auto p = data[row];
						w << p->x << p->y;
					}
					else {
						w << "" << "";
					}
				}

				w.Ln();
			}
		}
	};
}  // namespace HBGraph
