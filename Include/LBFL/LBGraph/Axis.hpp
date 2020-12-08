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

namespace LBGraph {
	/**
	 *	Axis objects map from work coordinates to screen coordinates.
	 */
	class Axis {
	protected:
		double min;		 //< The world-coordinate value corresponding to the lowest pixel value on the display axis covered by the plot area. For the x-axis, this will usually be the smaller of the interval endpoints covered by the plot area. For the y-axis, this will usually be the greater of the interval endpoints.

		double max;		 //< The world-coordinate value corresponding to the highest pixel value on the graph. For the x-axis, this will usually be the greater of the interval endpoints covered by the plot area. For the y-axis, this will usually be the lesser of the interval endpoints.

		bool isVisible;  //< Boolean indicator which is true if and only if the axis is to be rendered on the display.

	public:
		Axis(
			double min,			   //< The world-coordinate value corresponding to the lowest pixel value on the graph. For the x-axis, this will usually be the smaller of the interval endpoints covered by the plot area. For the y-axis, this will usually be the greater of the interval endpoints.

			double max,			   //< The world-coordinate value corresponding to the highest pixel value on the graph. For the x-axis, this will usually be the greater of the interval endpoints covered by the plot area. For the y-axis, this will usually be the lesser of the interval endpoints.

			bool isVisible = true  //< Boolean indicator which is true if and only if the axis is to be rendered on the display.
		) : min( min ), max( max ), isVisible( isVisible ) {}

		/**
		 * Get the world-coordinate value corresponding to the lowest pixel value on the display axis covered by the plot area.
		 * @returns this->min.
		 */
		double Min() {
			return min;
		}

		/**
		 * Set the world-coordinate value corresponding to the lowest pixel value on the display axis covered by the plot area.
		 * @param value The new value to be assigned to this->min.
		 */
		Axis &  SetMin( double value ) {
			min = value; return *this;
		}

		/**
		 * Get the world-coordinate value corresponding to the highest pixel value on the display axis covered by the plot area.
		 * @returns this->min.
		 */
		double Max() {
			return max;
		}

		/**
		 * Set the world-coordinate value corresponding to the highest pixel value on the display axis covered by the plot area.
		 * @param value The new value to be assigned to this->max.
		 */
		Axis &  SetMax( double value ) { max = value; return *this; }

		/**
		 * Get the current visibility of the axis.
		 *	@returns this->isVisible
		 */
		bool IsVisible() { return isVisible; }

		/**
		 * Sets the current visibility of the axis.
		 * @param value The new visibility indicator fo the axis (true == visible; false == invisible).
		 */
		Axis & SetVisible( bool value ) {
			isVisible = value;
			return *this;
		}

		/**
		 * Maps a world position to a screen pixel position.
		 * @param t Numeric value.
		 * @param lowPx The lower bound of the screen axis covered by the plot area (e.g. left or top).
		 * @param highPx The upper bound of the screen axis covered by the plot area (e.g. right or bottom).
		 * @returns The pixel location relative to subrange lowPx..highPx .
		 */
		virtual double ToScreen( double t, int lowPx, int highPx ) = 0;

		/**
		 * Maps a screen pixel position to a world position.
		 * @param px Pixel position.
		 * @param lowPx The lower bound of the screen axis covered by the plot area (e.g. left or top).
		 * @param highPx The upper bound of the screen axis covered by the plot area (e.g. right or bottom).
		 * @returns The numeric value that projects to px on the screen.
		 */

		virtual double ToWorld( double px, int lowPx, int highPx ) = 0;

		Axis & SetBounds( double lowPxValue, double highPxValue ) {
			this->min = lowPxValue;
			this->max = highPxValue;
			return *this;
		}
	};

	/**
	 * Specialised Axis which interpolates linearly between world and screen intervals.
	 */
	class LinearAxis : public Axis {
	public:
		/**
		 * Construct a linear axis.
		 */
		LinearAxis(
			double min,			   //< The world-coordinate value corresponding to the least pixel value in the screen interval.
			double max,			   //< The world-coordinate value corresponding to the greatest pixel value in the screen interval.
			bool isVisible = true  //< indicates whether the axis should be rendered.
		)
			: Axis( min, max, isVisible ) {}

		/**
		 * Maps a world position to a screen pixel position.
		 * @param t Numeric value.
		 * @param lowPx The lower bound of the screen axis covered by the plot area (e.g. left or top).
		 * @param highPx The upper bound of the screen axis covered by the plot area (e.g. right or bottom).
		 * @returns The pixel location relative to subrange lowPx..highPx .
		 */

		double ToScreen( double t, int lowPx, int highPx ) override {
			if ( lowPx > highPx ) {
				std::swap( lowPx, highPx );
			}
			return min == max ? NAN : lowPx + (t - min) * (highPx - lowPx) / (max - min);
		}

		/**
		 * Maps a screen pixel position to a world position.
		 * @param px Pixel position.
		 * @param lowPx The lower bound of the screen axis covered by the plot area (e.g. left or top).
		 * @param highPx The upper bound of the screen axis covered by the plot area (e.g. right or bottom).
		 * @returns The numeric value that projects to px on the screen.
		 */

		double ToWorld( double px, int lowPx, int highPx ) override {
			if ( lowPx > highPx ) {
				std::swap( lowPx, highPx );
			}
			return highPx == lowPx ? NAN : min + (px - lowPx) * (max - min) / (highPx - lowPx);
		}
	};

	/**
	 * Specialised Axis which interpolates logarithmically between world and screen intervals.
	 */
	class LogarithmicAxis : public Axis {
	public:
		/**
		 * Construct a logarithmic axis.
		 */
		LogarithmicAxis(
			double min,			   //< The world-coordinate value corresponding to the least pixel value in the screen interval.
			double max,			   //< The world-coordinate value corresponding to the greatest pixel value in the screen interval.
			bool isVisible = true  //< indivcates whether the axis should be rendered.
		)
			: Axis( min, max, isVisible ) {}

		/**
		 * Maps a world position to a screen pixel position.
		 * @param t Numeric value.
		 * @param lowPx The lower bound of the screen axis covered by the plot area (e.g. left or top).
		 * @param highPx The upper bound of the screen axis covered by the plot area (e.g. right or bottom).
		 * @returns The pixel location relative to subrange lowPx..highPx .
		 */

		double ToScreen( double t, int lowPx, int highPx ) override {
			if ( lowPx > highPx ) {
				std::swap( lowPx, highPx );
			}
			return max == min ? NAN : lowPx + (log( t ) - log( min )) * (highPx - lowPx) / (log( max ) - log( min ));
		}

		/**
		 * Maps a screen pixel position to a world position.
		 * @param px Pixel position.
		 * @param lowPx The lower bound of the screen axis covered by the plot area (e.g. left or top).
		 * @param highPx The upper bound of the screen axis covered by the plot area (e.g. right or bottom).
		 * @returns The numeric value that projects to px on the screen.
		 */

		double ToWorld( double px, int lowPx, int highPx ) override {
			if ( lowPx > highPx ) {
				std::swap( lowPx, highPx );
			}
			return highPx == lowPx ? NAN : exp( log( min ) + (px - lowPx) * (log( max ) - log( min )) / (highPx - lowPx) );
		}
	};
}  // namespace HBGraph
