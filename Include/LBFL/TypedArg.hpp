#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Group.H>
#include <string>
#include <vector>
#include "Requirement.hpp"
#include "PropertyChangedEventSource.hpp"

using namespace std;

namespace LBFL {
	/**
			<summary>
				A class of slightly reusable input controls for getting parameters such as numeric,
				string or file names.

				This control resembles Arg, however the internal implementation of the Value() function
				is slightly different.
				<list>
					<item>
						Arg uses a control which has a <c>char * value()</c>
					</item>
					<item>
						TypedArg uses a control which has a <c>T & value()</c>, so it is not necessary to
						parse the string out of the inner control.
					</item>
				</list>
			</summary>

			<typeparam name="T">
				The value data type.
			</typeparam>
			<typeparam name="FieldType">
				The data type for an input control that will be used to enter the value.
			</typeparam>
		 */
	template <typename T, typename FieldType>
	struct TypedArg : 
		public virtual Fl_Group, 
		public virtual PropertyChangedEventSource
	{
		Requirement isRequired;
		string name;
		string displayName;
		string help;
		Fl_Box leftBox;
		FieldType inputField;
		Fl_Button helpButton;
		const T defaultValue;

		TypedArg(const string &name, const T &defaultValue, const string &displayName, Requirement isRequired, const string &help,
			int labelWidth, int top, int w, int h
		) : Fl_Group(0, top, w, h),
			PropertyChangedEventSource(this),
			isRequired(isRequired),
			name(name),
			defaultValue( defaultValue ),
			displayName(displayName),
			help(help),
			leftBox(Fl_Boxtype::FL_FLAT_BOX, 0, top, w, h, ""),
			inputField(labelWidth + h + 5, top, w - labelWidth - h - 5, h),
			helpButton(labelWidth, top, h, h, "?")
			//
		{
			add(&leftBox);
			leftBox.label(this->displayName.c_str());
			leftBox.align(FL_ALIGN_LEFT | FL_ALIGN_INSIDE);
			add(&inputField);
			add(&helpButton);
			end();

			resizable(inputField);
			SetValue(defaultValue);

			helpButton.callback(HelpClicked, this);
			inputField.callback(InputChanged, this);
		}

		virtual ~TypedArg() {}

		/**
		 * Static event handler for FLTK change event on input field. Forwards a PropertyChanged message to each registered property changed event listener.
		 * @param sender the control which originated the change.
		 * @param context the host object containing the control.
		 */
		static void InputChanged(Control *sender, void *context)
		{
			// cerr << "TypedArg: InputChanged\n";
			auto arg = (TypedArg *)context;
			arg->NotifyPropertyChanged("Value");
		}

		void Reset() {
			SetValue( defaultValue );
		}

		/**
		 * Gets the current value of the embedded control.
		 * @returns The current value of the embedded control.
		 */
		T Value() const
		{
			auto val = inputField.value();
			return val;
		}

		/**
		 * Overwrites the current value of the embedded control.
		 * @param val a new value.
		 */
		void SetValue(const T & val)
		{
			if ( val == Value() ) return;

			inputField.value(val);
			NotifyPropertyChanged(NAMEOF(Value));
		}

		/**
		 * Gets the name of this control.
		 * @returns The name of this control.
		 */
		const string &Name()
		{
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
		static void HelpClicked(Control *sender, void *context)
		{
			fl_message("%s", ((TypedArg *)context)->help.c_str());
		}
	};
}  // namespace LBFL
