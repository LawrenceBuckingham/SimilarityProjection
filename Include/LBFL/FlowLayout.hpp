#pragma once

#include <FL/Fl.H>
#include <FL/fl_ask.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_Scroll.H>

#include <vector>

#include "Group.hpp"

namespace LBFL {

	/**
	 *	A group that implements something similar to java.awt.GridLayout.
	 *	Child controls are resized to become rectangles of equal width and height.
	 */
	class FlowLayout : public Group {
	protected:
		/** The horizontal gap (in pixels) to leave between columns. */
		int hGap;

		/** The vertical gap (measured in pixels) to leave between rows. */
		int vGap;

	public:
		/**
		 *	FlowLayout constructor.
		 *
		 *	<strong>NB:</strong>
		 *	<list>
		 *		<item>
		 *			Unlike pre-defined FLTK containers, this constructor automatically calls <c>end()</c>
		 *			to prevent it automatically snaffling newly created controls.
		 *		</item>
		 *		<item>
		 *			Use the <c>add()</c> method to insert child controls explicitly.
		 *		</item>
		 *	</list>
		 *</summary>
		 *<param name="x">
		 *	The nominal position of the left edge of the container.</param>
		 *<param name="y">
		 *	The nominal position of the top edge of the container.</param>
		 *<param name="w">
		 *	The nominal width of the container.</param>
		 *<param name="h">
		 *	The nominal height of the container.</param>
		 *<param name="hGap">
		 *	The horizontal gap (in pixels) to leave between columns.
		 *</param>
		 *<param name="vGap">
		 *	The vertical gap (in pixels) to leave between rows.</param>
		 */

		FlowLayout( int x, int y, int w, int h, int hGap = 5, int vGap = 5 ) : Group( x, y, w, h ), hGap( hGap ), vGap( vGap ) {
			end();
		}

		/**
		<summary>
			Resize the container, and lay out children to form a grid of equal&ndash;sized cells.
			<list>
				<item>
					This method hides <c>Fl_Group::add()</c>. Which hopefully will not be a problem.
				</item>
			</list>
		</summary>
		<param name="control">The address of the control to be added to this container.</param>
		 */

		void add( Fl_Widget * control ) {
			Fl_Group::add( control );
			ResizeImpl();
		}

		/**
		<summary>
			Resize the container, and lay out children to form a grid of equal&ndash;sized cells.
			<list>
				<item>
					This method hides <c>Fl_Group::add()</c>. Which hopefully will not be a problem.
				</item>
			</list>
		</summary>

		<param name="control">A reference to the control to be added to this container.</param>
		 */

		void add( Fl_Widget & control ) {
			Fl_Group::add( control );
			ResizeImpl();
		}

		/**
		<summary>
			Resize the container, and lay out children to form a grid of equal&ndash;sized cells.
			<para>
				<strong>NB:</strong> The real work is done by ResizeImpl.
			</para>
		</summary>

		<param name="xNew">The new position of the left edge of the container.</param>
		<param name="yNew">The new position of the top edge of the container.</param>
		<param name="wNew">The new width of the container.</param>
		<param name="hNew">The new height of the container.</param>
		*/

		void resize( int xNew, int yNew, int wNew, int hNew ) override {
			// cerr << "HERE 2!\n";
			Fl_Widget::resize( xNew, yNew, wNew, hNew );
			ResizeImpl();
		}

		/**
		<summary>
			Lay out children to form a grid of equal&ndash;sized cells.
		</summary>
		*/

		virtual void ResizeImpl() {
			//fprintf(stderr, "HERE 3: x = %d, y = %d, w = %d, h = %d, rows = %d, cols = %d\n", x(), y(), w(), h(), rows, cols);

			auto C = this->children();
			auto childArray = this->array();

			Fl_Boxtype b = box();
			int xx = x() + Fl::box_dx( b );  // was 9 instead of dx...
			int yy = y() + Fl::box_dy( b );
			int ww = w() - Fl::box_dw( b );
			int hh = h() - Fl::box_dh( b );

			int left = -hGap, top = -vGap, width = ww, height = hh;
			int currentLineHeight = 0;
			vector<vector<Fl_Widget *>> rows(1);
			vector<int> totalLineWidth(1);
			int lineCounter = 0;

			for ( int i = 0; i < C; i++ ) {
				auto child = childArray[i];

				if ( left + hGap + child->w() >= width ) {
					left = -hGap;
					top += vGap + currentLineHeight;
					currentLineHeight = 0;
					lineCounter++;
					totalLineWidth.resize( lineCounter + 1 );
					rows.resize( lineCounter + 1 );
				}

				if ( child->h() > currentLineHeight ) {
					currentLineHeight = child->h();
				}

				//fprintf(stderr, "Laying out child at (%d, %d, %d, %d)\n", left + x(), top + y(), width, height);

				child->resize( left + hGap + xx, top + vGap + yy, child->w(), child->h() );

				left += hGap + child->w();

				totalLineWidth[lineCounter] = left;
				rows[lineCounter].push_back( child );
			}

			for ( int l = 0; l < (int) totalLineWidth.size(); l++ ) {
				int adjustment =
					alignment == Align::Left ? 0 :
					alignment == Align::Centre ? (ww - totalLineWidth[l]) / 2 : ww - totalLineWidth[l];

				for ( auto child: rows[l] ) {
					child->position( child->x() + adjustment, child->y() );
				}
			}
		}

	public:
		enum Align { Left, Centre, Right };

	private:
		Align alignment = Align::Left;

	public:
		Align Alignment() {
			return alignment;
		}

		void SetAlignment( const Align val ) {
			if ( val == this->alignment ) return;

			this->alignment = val;
		}

	};
}
