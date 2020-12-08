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

	/**
		<summary>
			A class of slightly reusable input controls for getting parameters such as numeric,
			string or file names.
		</summary>

		<typeparam name="T">The value data type.</param>
		<typeparam name="FieldType">The data type for an input control that will be used to enter the value.</param>
		 */
	template <typename T, typename FieldType>
	class ValueWidget :
		public FieldType
		//
	{
	public:
		using ActionType = function<void()>;
	private:
		ActionType action;

		static ActionType DoNothing() {
			return []() {};
		}

	public:
		ValueWidget( int left, int top, int w, int h, const string & defaultValue, decltype(action) action = DoNothing() )
			: FieldType( left, top, w, h ), action(action)
			//
		{
			SetValue( defaultValue );
			Fl_Widget::callback( InputChanged, this );
		}

		virtual ~ValueWidget() {}

		static void InputChanged( Control *sender, void *context ) {
			auto arg = (ValueWidget *) context;
			arg->action();
		}

		/**
		 * Returns true if and only if the text field is non-empty. The value in the field may be garbage, but that is a different matter.
		 */
		bool HasValue() {
			string s = String::Trim( Fl_Input::value() );
			return s.length() > 0;
		}

		/**
		 * Gets the current value of the control.
		 *
		 * @returns The current value of the control.
		 * @throws FormatException if the text value is not valid for data type.
		 */
		T Value() const {
			return Util::Parse<T>( Fl_Input::value() );
		}

		/**
			 * Sets the value of the control.
			 *
			 * @param val The new value for the control.
			 */
		void SetValue( const T val ) {
			Fl_Input::value( Util::ToString<T>( val ).c_str() );
		}

		const decltype(action) & Action() const { return action; };

		void SetAction( const decltype(action) & value ) { action = value; }
	};
}
