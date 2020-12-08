#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Group.H>
#include <FL/fl_draw.H>
#include <memory>
#include <string>
#include <vector>
#include <String.hpp>
#include "PropertyChangedEventSource.hpp"
#include "Requirement.hpp"

using namespace QutBio;
using namespace std;

namespace LBFL {
#include "Util.hpp"
	/// Specialised grid layout which uniformly distributes a list of objects 
	///	having a required minimum size down the page.
	class VerticalFitLayout : public GridLayout {
	private:
		int minRowHeight;

	public:
		/**
		 *	Initialise a VerticalFitLayout object which will evenly distribute 
		 *	as many rows as possible, while ensuring that each row is at least 
		 *	minRowHeight pixels high. If minRowHeight is less than or equal to 
		 *	zero then a single row which dynamically occupies the full height 
		 *	of the container is maintained.
		 *
		 *	@param minRowHeight The minimum vertical height required for each 
		 *		row. To create a single row occupying all available space, use
		 *		0 or negative.
		 */
		VerticalFitLayout( int minRowHeight ) : minRowHeight( minRowHeight ) {}

		/**
		 *	Destructor.
		 */
		virtual ~VerticalFitLayout() {}

		/**
		 *	Override resize to recompute the required number of rows and resize 
		 *	child controls.
		 */
		void resize( int x, int y, int w, int h ) override {
			Fl_Boxtype b = box();
			int hh = h - Fl::box_dh( b );

			GridLayout::resize( x, y, w, h );
			int requiredRows = minRowHeight <= 0 ? 1 : std::max( 1, hh / minRowHeight );
			SetRows( requiredRows );
		}
	};
}
