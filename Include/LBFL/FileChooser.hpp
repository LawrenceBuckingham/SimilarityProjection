#pragma once

#include <FL/Fl.H>
#include <FL/fl_ask.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Double_Window.H>
#include <FL/Fl_Group.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Widget.H>
#include <FL/Fl_File_Chooser.H>

#include <vector>
#include <string>

#include "Exception.hpp"
#include "BorderLayout.hpp"
#include "PropertyChangedEventSource.hpp"

using namespace std;

namespace LBFL {

	/**
	 *	A group that implements a FileChooser which stashes its result in an (editable)
	 *	text area.
	 */

	class FileChooser : 
		public virtual BorderLayout, 
		public virtual PropertyChangedEventSource {
	public:
		enum class SelectionMode {
			/** Select a single, existing file. **/
			Single = Fl_File_Chooser::SINGLE,

			/** Select one or more existing files. **/
			Multi = Fl_File_Chooser::MULTI,

			/** Select a single existing file, or create a new one. **/
			Create = Fl_File_Chooser::CREATE,

			/** Select a single existing directory. */
			Directory = Fl_File_Chooser::DIRECTORY,
		};

	private:
		/** The text to display in the title bar. **/
		string title;

		/** A text field which shows the user which file was chosen. **/
		Fl_Input * display;

		/** A button that opens the proper file chooser **/
		Fl_Button * open;

		/** The file selection mode. **/
		SelectionMode mode;

		/** A copy of the value for use by customers. */
		string value_;

		string directory;
		string pattern;

	public:

		/**
			<summary>
				FileChooser constructor.

				<strong>Notes:</strong>
				<list>
					<item>
						The current Fl_Group is preserved intact;
					</item>
					<item>
						Multiple selection is not supported.
					</item>
				</list>
			</summary>
			<param name="left">The nominal position of the left edge of the control.</param>
			<param name="top">The nominal position of the top edge of the control.</param>
			<param name="width">The nominal width of the control.</param>
			<param name="height">The nominal height of the control.</param>
			<param name="mode">The selection mode. Note that SelectionMode::Multi is not supported.</param>
			<param name="title">The title bar text for the file chooser child window that opens.</param>
			<exception cref="QutBio::Exception">Thrown if multiple selection mode is requested.</exception>
		*/
		FileChooser(int left, int top, int width, int height,
			SelectionMode mode,
			const string & directory,
			const string & pattern,
			const string & title = "Please choose a file:"
		) :
			BorderLayout(left, top, width, height),
			PropertyChangedEventSource(this),
			directory(directory),
			pattern(pattern),
			title(title),
			mode(mode)
			//
		{
			if (mode == SelectionMode::Multi) {
				throw QutBio::Exception("Multiple selection is not supported.", FileAndLine);
			}

			end();

			Fl_Group::current(0);
			add(display = new Fl_Input(0, 0, 0, height), BorderLayout::Location::Centre);
			add(open = new Fl_Button(0, 0, 100, height, "Choose file..."), Location::East);

			ResizeImpl();
			
			open->callback(OpenClicked, this);
			display->callback(DisplayChanged, this);
		}

		/**
		<summary>
			Gets the currently selected file name.
		</summary>
		 */

		const string & value() const {
			return value_;
		}

		/**
		<summary>
			Sets the currently selected file name.

			This of course comes with no guarantee that the nominated file <em>exists</em>!
		</summary>
		<param name="val">A string which will become the initial selected file name.</param>
		 */

		void value(const string & val) {
			if (value_ != val) {
				value_ = val;
				this->display->value(val.c_str());
				// cerr << "FileChooser: Sending PropertyChanged: value = " << value_ << "\n";
				NotifyPropertyChanged(NAMEOF(value));
				auto cb = callback();
				if (cb) cb(this, this->user_data());
			}
		}

		/**
			<summary>
				Event handler for the open button. Displays a FLTK file chooser and stores the selected file (if any) in the containing FileChooser. It may not make any sense to call this, but just in case...
			</summary>
			<param name="sender">The control which is the immediate source of the message.</param>
			<param name="fileChooser">The address of a (assumed to be) FileChooser in which the selected file name will be saved.</param>
		 */

		static void OpenClicked(Fl_Widget * sender, void * fileChooser) {
			auto hbChooser = (FileChooser *)fileChooser;
			hbChooser->ShowDialog();
		}

		/// Callback for change event
		static void DisplayChanged(Fl_Widget * sender, void * fileChooser)
		{
			auto hbChooser = (FileChooser *)fileChooser;
			hbChooser->value(hbChooser->display->value());
		}

		/**
			<summary>
				Displays a FLTK file chooser in the designated mode, and stashes the
				resulting selected file name in this->value() if something was selected.
				If nothing was selected, then the value is preserved without change.
			</summary>
		 */
		void ShowDialog() {
			Fl_File_Chooser flChooser(directory.c_str(), pattern.c_str(), static_cast<int>(mode), title.c_str());
			flChooser.value(value().c_str());
			flChooser.show();

			while (flChooser.shown()) {
				Fl::wait();
			}

			if (flChooser.value()) {
				value(flChooser.value(1));
			}
		}

		const string & Pattern() const { return pattern; }
		void SetPattern( const string & val ) { pattern = val; }

		const string & Directory() const { return directory; }
		void SetDirectory( const string & val ) { directory = val; }
	};

	/**
		<summary>
			Specialisation of FileChooser which opens a single existing file.
		</summary>
	 */
	class InFileChooser : public virtual FileChooser {
	public : 
		/**
			<summary>
				InFileChooser constructor.

				<strong>Notes:</strong>
				<list>
					<item>
						The current Fl_Group is preserved intact;
					</item>
				</list>
			</summary>
			<param name="left">The nominal position of the left edge of the control.</param>
			<param name="top">The nominal position of the top edge of the control.</param>
			<param name="width">The nominal width of the control.</param>
			<param name="height">The nominal height of the control.</param>
			<param name="title">The title bar text for the file chooser child window that opens.</param>
		 */
		InFileChooser(int left, int top, int width, int height,
			const string & directory = ".",
			const string & pattern = "*",
			const string & title = "Please choose a file:"
		) : BorderLayout(left, top, width, height),
			PropertyChangedEventSource(this),
			FileChooser( left, top, width, height, FileChooser::SelectionMode::Single, directory, pattern, title) {}
	};

	/**
			Specialisation of FileChooser which opens a single existing file, or 
			allows a new one to be created.
	 */
	class OutFileChooser : public virtual FileChooser {
	public : 
		/**
			<summary>
				OutFileChooser constructor.

				<strong>Notes:</strong>
				<list>
					<item>
						The current Fl_Group is preserved intact;
					</item>
				</list>
			</summary>
			<param name="left">The nominal position of the left edge of the control.</param>
			<param name="top">The nominal position of the top edge of the control.</param>
			<param name="width">The nominal width of the control.</param>
			<param name="height">The nominal height of the control.</param>
			<param name="title">The title bar text for the file chooser child window that opens.</param>
		 */
		OutFileChooser(int left, int top, int width, int height,
			const string & directory = ".",
			const string & pattern = "*",
			const string & title = "Please choose a file:"
		) : BorderLayout(left, top, width, height),
			PropertyChangedEventSource(this),
			FileChooser( left, top, width, height, FileChooser::SelectionMode::Create, directory, pattern, title) {}
	};
}
