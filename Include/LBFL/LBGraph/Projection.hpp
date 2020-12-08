#pragma once

#include <functional>

namespace LBGraph {
/**
 *	Delegate type which projects from world coordinates to screen coordinates, returning true iff the point is within bounds of the current viewport.
 */
using Projection = function<bool(double wx, double wy, int &sx, int &sy)>;
}