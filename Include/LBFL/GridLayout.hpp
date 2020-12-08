#pragma once

#include <FL/Fl.H>
#include <FL/fl_ask.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_Scroll.H>
#include "Group.hpp"

#include <vector>

namespace LBFL {

	/**
	 *	A group that implements something similar to java.awt.GridLayout.
	 *	Child controls are resized to become rectangles of equal width and height.
	 */
	class GridLayout : public Group {
	protected:
		/** The number of rows in the grid. */
		int rows;

		/** The number of columns in the grid. */
		int cols;

		/** The horizontal gap (in pixels) to leave between columns. */
		int hGap;

		/** The vertical gap (measured in pixels) to leave between rows. */
		int vGap;

	public:
		/**
		<summary>
			GridLayout constructor.

			<strong>NB:</strong>
			<list>
				<item>
					Unlike pre-defined FLTK containers, this constructor automatically calls <c>end()</c>
					to prevent it snaffling newly created controls.
				</item>
				<item>
					Use the <c>add()</c> method to insert child controls explicitly.
				</item>
			</list>
		</summary>
		<param name="x">
			The nominal position of the left edge of the container.</param>
		<param name="y">
			The nominal position of the top edge of the container.</param>
		<param name="w">
			The nominal width of the container.</param>
		<param name="h">
			The nominal height of the container.</param>
		<param name="rows">
			The number of rows.</param>
		<param name="cols">
			The number of columns.</param>
		<param name="hGap">
			The horizontal gap (in pixels) to leave between columns.
		</param>
		<param name="vGap">
			The vertical gap (in pixels) to leave between rows.</param>
		*/
		
		GridLayout(
			int x = 0, 
			int y = 0, 
			int w = 0, 
			int h = 0, 
			int rows = 1, 
			int cols = 1, 
			int hGap = 5, 
			int vGap = 5
		) : Group(x, y, w, h), rows(rows), cols(cols) {
			// cerr << "HERE 1!\n";
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

		void add(Fl_Widget * control) {
			Fl_Group::add(control);
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

		void add(Fl_Widget & control) {
			Fl_Group::add(control);
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

		void resize(int xNew, int yNew, int wNew, int hNew) override {
			// cerr << "HERE 2!\n";
			Fl_Widget::resize(xNew, yNew, wNew, hNew);
			ResizeImpl();
		}

		/**
		<summary>
			Lay out children to form a grid of equal&ndash;sized cells.
		</summary>
		*/

		virtual void ResizeImpl() {
			//fprintf(stderr, "HERE 3: x = %d, y = %d, w = %d, h = %d, rows = %d, cols = %d\n", x(), y(), w(), h(), rows, cols);

			Fl_Boxtype b = box();
			int xx = x() + Fl::box_dx( b );  // was 9 instead of dx...
			int yy = y() + Fl::box_dy( b );
			int ww = w() - Fl::box_dw( b );
			int hh = h() - Fl::box_dh( b );

			auto C = this->children();
			auto childArray = this->array();
			auto requiredRows = std::max(rows, (C + cols - 1) / cols);

			for (int i = 0; i < C; i++) {
				auto col = i % cols;
				auto right = ww * (col + 1) / cols;
				auto left = ww * col / cols;
				auto width = right - left;

				auto row = i / cols;
				auto top = hh * row / rows;
				auto bottom = hh * (row + 1) / rows;
				auto height = bottom - top;

				//fprintf(stderr, "Laying out child at (%d, %d, %d, %d)\n", left + x(), top + y(), width, height);

				childArray[i]->resize(left + xx, top + yy, width, height);
			}
		}

		int Rows() { return rows; }

		void SetRows( int value ) {
			if ( this->rows == value ) return;

			this->rows = value;
			ResizeImpl();
		}

		int Cols() { return rows; }

		void SetCols( int value ) {
			if ( this->cols == value ) return;

			this->cols = value;
			ResizeImpl();
		}

		void SetRowsCols( int rowValue, int colValue ) {
			if ( rows == rowValue && cols == colValue ) return;

			rows = rowValue;
			cols = colValue;
			ResizeImpl();
		}
	};
}
