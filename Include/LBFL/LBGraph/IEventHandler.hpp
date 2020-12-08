#pragma once

#include <FL/Fl_Widget.H>

namespace LBGraph {
	class IEventHandler {
		 public:
		 virtual ~IEventHandler() {}

		 virtual bool HandleEvent( Fl_Widget * src, int eventCode ) = 0;
	 };
 }