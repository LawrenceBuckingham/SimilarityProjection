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
using namespace std;

namespace LBFL {

	/**
	 *	A group that implements something similar to java.awt.BorderLayout.
	 */
	class BorderLayout : public Group {
	protected:
		/** Track the controls resident in each location (independently of the underlying FLTK child collection, which does not suit my purposes well.) */
		std::vector<std::vector<Fl_Widget *>> children;

	public:
		/**
		 *	The 5 locations know to BorderLayout.
		 */
		enum class Location {
			/** Stretch to fill the central unoccupied region. */
			Centre,

			/** Stretch horizontally across the top of the container, using the existing height of the control. */
			North,

			/** Stretch horizontally across the bottom of the container, using the existing height of the control. */
			South,

			/** Stretch vertically up the right-hand edge of the container, using the existing width of the control. */
			East,

			/** Stretch vertically up the left-hand edge of the container, using the existing width of the control. */
			West
		};

		/**
		 *	<summary>
		 *	Constructs a new BorderLayout, which emulates a Java Container managed
		 *	by a java.awt.BorderLayout layout manager.
		 *
		 *	<b>NB:</b> Unlike the FLTK groups, this container will not automatically
		 *	snaffle newly constructed Widgets; the reason is that we cannot reliably
		 *	<code>add()</code> widgets without specifying the location parameter.
		 *	</summary>
		 *
		 *	<param name="x">The left edge of the control, nominally relative to the window.</param>
		 *	<param name="y">The top edge of the control, nominally relative to the window.</param>
		 *	<param name="w">The nominal width of the control.</param>
		 *	<param name="h">The nominal height of the control.</param>
		 */
		BorderLayout(int x = 0, int y = 0, int w = 0, int h = 0) : Group(x, y, w, h), children(5) {
			end();
		}

		/**
		 *	<summary>Adds a child control at the centre of the container.
		 *	<list>
		 *	<item>This hides the <tt>Fl_Group::add</tt> method.</item>
		 *	</list>
		 *	</summary>
		 *	<param name="control">Reference to the new child control to add.</param>
		 */
		void add(Fl_Widget & control) /*override*/ {
			add(&control, Location::Centre);
		}

		/**
		 *	<summary>Adds a child control at the centre of the container.
		 *	<list>
		 *	<item>This hides the <tt>Fl_Group::add</tt> method.</item>
		 *	</list>
		 *	</summary>
		 *	<param name="control">Address of the new child control to add.</param>
		 */

		void add(Fl_Widget * control) /*override*/ {
			add(control, Location::Centre);
		}

		/**
		 *	<summary>
		 *	Adds a child control at a nominated location in the container.
		 *	</summary>
		 *	<param name="control">Address of the new child control to add.</param>
		 *	<param name="location">The Location where the child should be displayed.</param>
		 */

		virtual void add(Fl_Widget * control, Location location) {
			Fl_Group::add(control);
			auto & childList = children[static_cast<int>(location)];
			childList.insert(childList.begin(), control);
			ResizeImpl();
			redraw();
		}

		/**
			<summary>
				Defines new location and dimensions for the container, and lays out the child
				controls using nominal dimensions obtained by querying their <c>w()</c> and
				<c>h()</c> properties.
				<list>
					<item>Overrides fl_Group::resize.</item>
					<item>
						Uses <c>ResizeImpl</c> to reposition and size each child based on its
						location and innate dimensions.
					</item>
				</list>
			</summary>
			<param name="xNew">The new left edge of the container.</param>
			<param name="yNew">The new top edge of the container.</param>
			<param name="wNew">The new width of the container.</param>
			<param name="hNew">The new height of the container.</param>
		 */

		void resize(int xNew, int yNew, int wNew, int hNew) override {
			Fl_Widget::resize(xNew, yNew, wNew, hNew); // make new xywh values visible for children
			ResizeImpl();
		}

		/** Gets the address of the control in the North location. */
		std::vector<Fl_Widget *> & North() { 
			return children[static_cast<int>(Location::North)]; 
		}

		/** Gets the address of the control in the North location. */
		std::vector<Fl_Widget *> & South() {
			return children[static_cast<int>(Location::South)];
		}

		/** Gets the address of the control in the North location. */
		std::vector<Fl_Widget *> & East() {
			return children[static_cast<int>(Location::East)];
		}

		/** Gets the address of the control in the North location. */
		std::vector<Fl_Widget *> & West() {
			return children[static_cast<int>(Location::West)];
		}

		/** Gets the address of the control in the North location. */
		std::vector<Fl_Widget *> & Centre() {
			return children[static_cast<int>(Location::Centre)];
		}

		/**
			<summary>
				Lays out the child controls using nominal dimensions obtained by querying their <c>w()</c> and
				<c>h()</c> properties.
			</summary>
		 */
		virtual void ResizeImpl() {
			int cx, cy, cw, ch;

			Fl_Boxtype b = box();
			int xx = x() + Fl::box_dx( b );  // was 9 instead of dx...
			int yy = y() + Fl::box_dy( b );
			int ww = w() - Fl::box_dw( b );
			int hh = h() - Fl::box_dh( b );

			auto north = North();
			auto south = South();
			auto east = East();
			auto west = West();
			auto centre = Centre();

			auto northHeight = north.size() > 0 ? north[0]->h() : 0;
			auto southHeight = south.size() > 0 ? south[0]->h() : 0;
			auto westWidth = west.size() > 0 ? west[0]->w() : 0;
			auto eastWidth = east.size() > 0 ? east[0]->w() : 0;

			if (west.size() > 0) {
				cx = xx;
				cy = yy + northHeight;
				cw = westWidth;
				ch = hh - southHeight - northHeight;
				for ( auto child: west) {
					child->resize(cx, cy, cw, ch);
				}
			}

			if (east.size() > 0) {
				cx = xx + ww - eastWidth;
				cy = yy + northHeight;
				cw = eastWidth;
				ch = hh - southHeight - northHeight;
				for ( auto child: east) {
					child->resize(cx, cy, cw, ch);
				}
			}

			if (north.size() > 0) {
				cx = xx;
				cy = yy;
				cw = ww;
				ch = northHeight;
				for ( auto child: north) {
					child->resize(cx, cy, cw, ch);
				}
			}

			if (south.size() > 0) {
				cx = xx;
				cy = yy + hh - southHeight;
				cw = ww;
				ch = southHeight;
				for ( auto child: south) {
					child->resize(cx, cy, cw, ch);
				}
			}

			if (centre.size() > 0) {
				cx = xx + westWidth;
				cy = yy + northHeight;
				cw = ww - westWidth - eastWidth;
				ch = hh - northHeight - southHeight;
				for ( auto child: centre) {
					child->resize(cx, cy, cw, ch);
				}
			}
		}

		void remove( Fl_Widget & child ) {
			remove(&child);
		}

		void remove( Fl_Widget * child ) {
			Fl_Group::remove( *child );

			for ( auto & list: children ) {
				auto iter = list.begin();

				while ( iter != list.end() ) {
					if ( *iter == child ) {
						list.erase(iter);
						break;
					}
				}
			}
		}
	};
}
