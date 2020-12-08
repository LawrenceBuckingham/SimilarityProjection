#pragma once

#include <FL/Fl.H>
#include <Fl/Fl_Button.H>
#include <FL/Fl_Int_Input.H>
#include <FL/fl_draw.H>

#include <FlowLayout.hpp>
#include <BorderLayout.hpp>
#include <GridLayout.hpp>
#include <SparseSignature.hpp>
#include <ScrollArea.hpp>
#include <ValueWidget.hpp>

#include <FastaSequence.hpp>
#include <SparseSignature.hpp>
#include <ActionWidget.hpp>
#include <Box.hpp>
#include <ScatterPlot.hpp>
#include <VerticalFitLayout.hpp>
#include "Helpers.hpp"
#include "SigRank.hpp"
#include "DefaultChartHandler.hpp"

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;
using namespace LBGraph;

class PrecisionRecallSummary : public BorderLayout {
private:
	Box labelBox;
	string heading;
	ScatterPlot chart;
	LineSpec precisionLine, boundLine, meanLine;
	Series precision, lowerBound, upperBound, meanSeries;
	DefaultChartHandler handler;

	const PrecisionRecallStats &stats;
public:
	virtual ~PrecisionRecallSummary() {}

	PrecisionRecallSummary( decltype(stats) stats ) :
		stats( stats ),
		labelBox( 0, 0, 30, 30 ),
		precisionLine( FL_GREEN, 2, FL_SOLID ),
		precision( 0, &precisionLine, 0, "precision" ),
		boundLine( FL_GRAY, 2, FL_DASH ),
		lowerBound( 0, &boundLine, 0, "lower bound" ),
		upperBound( 0, &boundLine, 0, "upper bound" ),
		meanLine( FL_BLUE, 2, FL_DASH ),
		meanSeries( 0, &meanLine, 0, "mean precision"),
		handler(&chart)
		//
	{
		add( &labelBox, Location::North );
		add( &chart, Location::Centre );

		chart.Add( &lowerBound )
			.Add( &upperBound )
			.Add( &precision )
			.Add( &meanSeries )
			.SetMargin( 40, 10, 30, 30 )
			.SetAxesCross( 0, 0 )
			.SetAxisBounds( 0, 1, 1, 0 )
			.SetTickMarks( 10, "%0.1f", 10, "%0.1f" )
			.SetFillColour( FL_WHITE )
			.AddMouseHandler(&handler);

		labelBox.SetFontFamily( FL_COURIER )
			.SetFontSize( 14 )
			.Anchor( 0, 0.5, 5, 0, 0, 0.5 )
			.SetFillColour( FL_WHITE )
			.SetTextColour( FL_BLACK );
	}

	void GainFocus() {
		Refresh();
	}

	void Clear() {
		precision.Clear();
		lowerBound.Clear();
		upperBound.Clear();
		meanSeries.Clear();
		heading = "Query dataset size: ?, MAP: ?";
		labelBox.SetLabel( heading );
		redraw();
	}

	void Refresh() {
		precision.Clear();
		lowerBound.Clear();
		upperBound.Clear();
		meanSeries.Clear();

		const uint N = stats.prec.size();

		if ( N > 0 ) {
			vector<double> column( N );

			for ( uint i = 0; i <= stats.numSteps; i++ ) {
				double recall = double( i ) / stats.numSteps;

				for ( uint j = 0; j < N; j++ ) {
					column[j] = stats.prec[j][i];
				}

				std::sort( column.begin(), column.end() );

				double median = (column[N / 2] + column[(N - 1) / 2]) / 2;

				auto idx25 = std::min( std::max( (N * 25 + 50) / 100 - 1, uint(0) ), N-1 );
				double low = column[idx25];

				auto idx75 = std::min( std::max( (N * 75 + 50) / 100 - 1, uint(0) ), N-1 );
				double high = column[idx75];

				auto sum = std::accumulate( column.begin(), column.end(), 0.0 );

				precision.Add( recall, median );
				lowerBound.Add( recall, low );
				upperBound.Add( recall, high );
				meanSeries.Add( recall, sum / N );
			}
		}

		heading = String::Format( "Query dataset size: %d, MAP: %0.2f", N, stats.meanAveragePrecision );
		labelBox.SetLabel( heading );

		redraw();
	}
};
