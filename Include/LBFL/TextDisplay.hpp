#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Text_Display.H>

#include <cstdarg>
#include <cstdio>

namespace LBFL {

	/**
	 *	Enhances Fl_Text_Display by adding formatted text insertion like printf.
	 */

	class TextDisplay : public Fl_Text_Display {
		Fl_Text_Buffer buff;

	public:
		/**
		 *	TextDisplay constructor.
		 *
		 *	@param[in] left		The nominal position of the left edge of the control.
		 *	@param[in] top		The nominal position of the top edge of the control.
		 *	@param[in] width	The nominal width of the control.
		 *	@param[in] height	The nominal height of the control.
		 *	@param[in] label	The title bar text for the file chooser child window that opens.
		 */

		TextDisplay(
			int left,
			int top,
			int width,
			int height,
			const char * label = 0
		) :
			Fl_Text_Display( left, top, width, height, label ) {
			end();
			buffer( &buff );
		}

		~TextDisplay() {
			buffer(0);
		}

		/**
		 *	Inserts formatted data at the current insertion point.
		 *	@param[in]	format	A standard C format string.
		 *	@param[in]	...		Additional arguments as required to satisfy the format specifier.
		 *	@returns			On success, returns the number of characters written.
		 *						Otherwise, returns a negative value.
		 */

		int operator()( const char * format, ... ) {
			va_list args;
			va_start( args, format );
			auto result = operator()( format, args );
			va_end( args );
			return result;
		}

		/**
		 *	Inserts formatted data at the current insertion point.
		 *	@param[in]	format	A standard C format string.
		 *	@param[in]	args	Additional arguments as required to satisfy the format specifier.
		 *	@returns			On success, returns the number of characters written.
		 *						Otherwise, returns a negative value.
		 */

		int operator()( const char * format, va_list args ) {
			char * buffer = 0;
			size_t size = 0;
			FILE * f = open_memstream( &buffer, &size );

			if ( f == 0 ) {
				return -1;
			}

			int charsWritten = vfprintf( f, format, args );
			fclose( f );
			insert( buffer );
			free( buffer );

			return charsWritten;
		}

		/**
		 * Removes all characters fromt he text display.
		 */
		void Clear() {
			mBuffer->remove( 0, mBuffer->length() );
		}

		/**
		 * Hide the original member, which deletes all child controls from the object.
		 */
		void clear() {
			Clear();
		}
	};
}
