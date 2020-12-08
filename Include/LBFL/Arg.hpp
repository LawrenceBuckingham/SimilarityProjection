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
	class Arg :
		public virtual Fl_Group,
		public virtual PropertyChangedEventSource {
	public:
		Requirement isRequired;
		string name;
		string displayName;
		string help;
		Fl_Box leftBox;
		FieldType inputField;
		Fl_Button helpButton;
		const T defaultValue;

		Arg( const string &name, const T &defaultValue, const string &displayName, Requirement isRequired, const string &help,
			int labelWidth, int top, int w, int h
		) : Fl_Group( 0, top, w, h ),
			PropertyChangedEventSource( this ),
			isRequired( isRequired ),
			name( name ),
			defaultValue(defaultValue),
			displayName( displayName ),
			help( help ),
			leftBox( Fl_Boxtype::FL_FLAT_BOX, 0, top, w, h, "" ),
			inputField( labelWidth + h + 5, top, w - labelWidth - h - 5, h ),
			helpButton( labelWidth, top, h, h, "?" )
			//
		{
			add( &leftBox );
			leftBox.label( this->displayName.c_str() );
			leftBox.align( FL_ALIGN_LEFT | FL_ALIGN_INSIDE );
			add( &inputField );
			add( &helpButton );
			end();

			resizable( inputField );
			SetValue( defaultValue );

			helpButton.callback( HelpClicked, this );
			inputField.callback( InputChanged, this );
		}

		virtual ~Arg() {}

		void Reset() {
			SetValue(defaultValue);
		}

		static void InputChanged( Control *sender, void *context ) {
			auto arg = (Arg *) context;
			arg->NotifyPropertyChanged( "Value" );
		}

		/**
		 * Returns true if and only if the text field is non-empty. The value in the field may be garbage, but that is a different matter.
		 */
		bool HasValue() {
			string s = String::Trim( inputField.value() );
			return s.length() > 0;
		}

		/**
		 * Gets the current value of the control.
		 *
		 * @returns The current value of the control.
		 * @throws FormatException if the text value is not valid for data type.
		 */
		T Value() const {
			return Util::Parse<T>( inputField.value() );
		}

		/**
			 * Sets the value of the control.
			 *
			 * @param val The new value for the control.
			 */
		void SetValue( const T val ) {
			inputField.value( Util::ToString<T>( val ).c_str() );
		}

		/**
			 * Get the name of the control.
			 * @returns A reference to the string containing the name of the control.
			 */
		const string &Name() {
			return name;
		}

		/**
		<summary>
			Callback function which handles click events on the "help" button.

			The function displays the help message associated with the context
			object.
		</summary>
		<param name="sender">
			The control that instigated the call-back. This will probably be the help button.
		</param>
		<param name="context">
			The address of the context object: an instance of Arg which will
			produce the desired help string.
		</param>
		*/
		static void HelpClicked( Control *sender, void *context ) {
			fl_message( "%s", ((Arg *) context)->help.c_str() );
		}
	};
}  // namespace LBFL
