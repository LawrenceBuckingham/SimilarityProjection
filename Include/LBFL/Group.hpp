#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Group.H>

namespace LBFL {
	class Group : public Fl_Group {
	public:
		Group( int x = 0, int y = 0, int w = 0, int h = 0, const char * label = 0 ) :
			Fl_Group( x, y, w, h, label )
			//
		{
			end();
		}

		virtual ~Group() {
			// Remove all child controls to prevent the Fl_Group destructor from DELETING them all.
			// Use proper RAII to manage child controls.
			clear();
		}

		void clear() {
			int n;
			while ( (n = children()) > 0 ) {
				remove( child( n - 1 ) );
			}
			Fl_Group::clear();
		}

		void GetAvailableSize( int & width, int & height ) {
			Fl_Boxtype b = box();
			width = w() - Fl::box_dw( b );
			height = h() - Fl::box_dh( b );
		}
	};
}
