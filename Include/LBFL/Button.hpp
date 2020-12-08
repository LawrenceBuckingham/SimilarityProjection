#pragma once

#include <functional>
#include <FL/Fl.H>
#include <FL/Fl_Button.H>

using namespace std;

namespace LBFL {
	class Button : public Fl_Button {
	public:
		using ActionType = function<void( Button * )>;
	private:
		ActionType action;
		void * data;

		static ActionType DoNothing() {
			return [](Button *b) {};
		}

	public:
		Button(
			int x = 0,
			int y = 0,
			int w = 0,
			int h = 0,
			const char * label = 0,
			ActionType action = DoNothing(),
			void *data = 0
		) :
			Fl_Button( x, y, w, h, label ),
			action( action ),
			data( data )
			//
		{
			callback( Callback );
		}

		static void Callback( Fl_Widget * widget ) {
			auto button = (Button *) widget;
			button->action( button );
		}

		ActionType Action() const { return action; }

		void SetAction( ActionType value ) {
			action = value;
		}
	};
}
