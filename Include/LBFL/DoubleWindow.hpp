#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Double_Window.H>

namespace LBFL {
	struct DoubleWindow :
		public Fl_Double_Window
		// 
	{
		DoubleWindow( int w, int h, const char * label = 0 )
			: Fl_Double_Window( w, h, label )
			//
		{
			end();
		}

		virtual ~DoubleWindow() {
			clear();
		}

		void clear() {
			int n;
			while ( (n = children()) > 0 ) {
				remove( child( n - 1 ) );
			}
		}
	};
}
