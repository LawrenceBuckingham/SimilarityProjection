#pragma once

#include <FL/Fl.H>
#include <FL/Enumerations.H>
#include <FL/Fl_File_Chooser.H>
#include <ScatterPlot.hpp>
#include <cctype>

using namespace LBFL;

struct DefaultChartHandler : public LBGraph::ScatterPlot::MouseHandler {
	ScatterPlot * chart;

	DefaultChartHandler( ScatterPlot * chart ) : chart( chart ) {}

	virtual ~DefaultChartHandler() {}

	virtual bool Handle( const LBGraph::ScatterPlot::MouseEvent & event ) override {
		// fprintf( stderr, "%s(%d) %d\n", FileAndLine, event.eventCode );

		if ( event.source != chart ) return false;

		if ( event.eventCode == FL_KEYUP ) {
			auto eventKey = Fl::event_key();
			auto eventState = Fl::event_state();

			// fprintf( stderr, "%s(%d) (%c)(%d)\n", FileAndLine, eventKey, eventState );

			if ( tolower( eventKey ) == 's' ) {
				// fprintf( stderr, "%s(%d) (%c)(%d)\n", FileAndLine, eventKey, eventState );

				if ( eventState & (FL_CONTROL | FL_COMMAND) ) {
					auto outFile = ShowDialog();

					if ( outFile.length() > 0 ) {
						event.source->SaveToFile( outFile );
					}

					return true;
				}
			}
		}

		return false;
	}

	static string ShowDialog() {
		Fl_File_Chooser flChooser( ".", "*.chart.csv", Fl_File_Chooser::CREATE, "Save chart to file." );
		flChooser.show();

		while ( flChooser.shown() ) {
			Fl::wait();
		}

		if ( flChooser.value() ) {
			return flChooser.value( 1 );
		}
		else {
			return "";
		}
	}

};
