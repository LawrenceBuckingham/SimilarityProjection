#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Scroll.H>

namespace LBFL {
	class ScrollArea : public Fl_Scroll {
	public:
		ScrollArea( int x, int y, int w, int h, const char * label = 0 ) :
			Fl_Scroll( x, y, w, h, label )
			//
		{
			end();
		}

		virtual ~ScrollArea() {
			// Remove all child controls to prevent the Fl_Group destructor from DELETING them all.
			// Use proper RAII to manage child controls.
			clear();
		}

		void clear() {
			int n;
			while ( (n = children()) > 0 ) {
				remove( child( n - 1 ) );
			}
		}

		void resize ( int x, int y, int w, int h ) override {
			// fprintf(stderr, "ScrollArea resizing to %d, %d, %d, %d\n", x, y, w, h);
			Fl_Scroll::resize(x, y, w, h);
		}
	};
}
