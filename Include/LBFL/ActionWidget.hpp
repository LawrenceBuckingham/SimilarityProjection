#pragma once

#include <functional>
#include <FL/Fl.H>
#include <FL/Fl_Int_Input.H>

using namespace std;

namespace LBFL {
	template<typename Widget>
	class ActionWidget : public Widget {
	public:
		using ActionType = function<void()>;
	private:
		ActionType action;

		static ActionType DoNothing() {
			return []() {};
		}

	public:
		ActionWidget(
			int x = 0,
			int y = 0,
			int w = 0,
			int h = 0,
			const char * label = 0,
			ActionType action = DoNothing()
		) :
			Widget( x, y, w, h, label ),
			action( action )
			//
		{
			Fl_Widget::callback( Callback );
		}

		static void Callback( Fl_Widget * widget ) {
			auto control = (ActionWidget *) widget;
			control->action();
		}

		ActionType Action() const { return action; }

		void SetAction( ActionType value ) {
			action = value;
		}
	};
}
