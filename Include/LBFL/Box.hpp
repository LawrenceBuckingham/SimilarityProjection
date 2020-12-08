#pragma once

#include <functional>
#include <FL/Fl.H>
#include <FL/Fl_Box.H>
#include <FL/Enumerations.H>
#include <FL/Fl_Widget.H>
#include <FL/fl_draw.H>

namespace LBFL {
	/// A box with a label.
	class Box : public Fl_Box {
		string label_;

	public:
		Box(
			int x = 0,
			int y = 0,
			int w = 0,
			int h = 0,
			const string & label = ""
		) :
			Fl_Box( x, y, w, h, 0 ),
			label_(label)
			//
		{}

		Box( string & label ) : Fl_Box( 0, 0, 0, 0, 0 ), label_(label) {
			AutoFit( 0, 0, 0, 0 );
		}

		void draw() override {
			draw_box();
			Fl_Boxtype b = box();
			int xx = x() + Fl::box_dx( b );  // was 9 instead of dx...
			int yy = y() + Fl::box_dy( b );
			int ww = w() - Fl::box_dw( b );
			int hh = h() - Fl::box_dh( b );
			fl_push_clip( xx, yy, ww, hh );
			fl_color( color() );
			fl_rectf( xx, yy, ww, hh );
			DrawLabel();
			fl_pop_clip();
		}

		/// Horizontal displacement (in pixels) of label from computed anchor point
		int offsetX = 0;

		/// Vertical displacement (in pixels) of label from computed anchor point 
		int offsetY = 0;

		/// Relative horizontal position of anchor point within text bounding box.
		double anchorX = 0.5;

		/// Relative vertical position of anchor point within text bounding box.
		double anchorY = 0.5;

		/// Relative horizontal position of anchor point within this Box.
		double anchorXto = 0.5;

		/// Relative vertical position of anchor point within this Box.
		double anchorYto = 0.5;

		/// Draws the label, using anchor points rather than the FLTK mechanism, giving me improved control.
		void DrawLabel() {
			auto lab = label_.c_str();

			int tw = 0, th = 0;
			fl_font( labelfont(), labelsize() );
			fl_measure( lab, tw, th );

			int labelX = x() + w() * anchorXto + offsetX - anchorX * tw;
			int labelY = y() + h() * anchorYto + offsetY - anchorY * th - fl_descent() + fl_height();
			fl_color( labelcolor() );
			fl_draw( lab, labelX, labelY );
		}

		/**
		 *	Define the anchor points for the label within the containing box.
		 *
		 *	@param anchorX Relative horizontal position of anchor point within text bounding box.
		 *	@param anchorY Relative vertical position of anchor point within text bounding box.
		 *	@param offsetX Horizontal displacement (in pixels) of label from computed anchor point
		 *	@param offsetY Vertical displacement (in pixels) of label from computed anchor point
		 *	@param anchorXTo Relative horizontal position of anchor point within this Box.
		 *	@param anchorYTo Relative vertical position of anchor point within this Box.
		 */
		Box & Anchor( double anchorX, double anchorY, int offsetX, int offsetY, double anchorXTo, double anchorYTo ) {
			this->anchorX = anchorX;
			this->anchorY = anchorY;
			this->offsetX = offsetX;
			this->offsetY = offsetY;
			this->anchorXto = anchorXTo;
			this->anchorYto = anchorYTo;
			return *this;
		}

		/**
		 *	Gets the dimensions of the label.
		 */
		Pair<int, int> GetTextSize() const {
			Pair<int, int> res{ 0,0 };
			fl_font( labelfont(), labelsize() );
			fl_measure( label_.c_str(), res.item1, res.item2 );
			return res;
		}

		/**
		 *	Resizes and re-anchors the text so that it fits in the box with the designated margins. 
		 */
		Box & AutoFit( int marginLeft, int marginTop, int marginRight, int marginBottom ) {
			auto dim = GetTextSize();
			Fl_Boxtype b = box();
			w( dim.item1 + marginLeft + marginRight + Fl::box_dw( b ) );
			h( dim.item2 + marginTop + marginBottom + Fl::box_dh( b ) );
			Anchor( 0, 0, marginLeft, marginTop, 0, 0 );
			return *this;
		}

		Box & SetLabel( const string & s ) {
			this->label_ = s;
			return *this;
		}

		Box & SetFillColour( Fl_Color fillColour ) {
			this->color(fillColour);
			return *this;
		}

		Box & SetFontFamily( Fl_Font fontNumber ) {
			labelfont( fontNumber );
			return *this;
		}

		Box & SetFontSize( int fontSize ) {
			labelsize( fontSize );
			return *this;
		}

		Box & SetTextColour( Fl_Color textColour ) {
			labelcolor(textColour);
			return *this;
		}

		void label( const char * s ) {
			throw NotImplementedException(FileAndLine);
		}
	};
}
