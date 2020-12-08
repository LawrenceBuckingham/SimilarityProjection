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
#include "DefaultChartHandler.hpp"

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;
using namespace LBGraph;

class PrecisionRecallChart : public ListItem<::Ranking> {
private:
	Box labelBox;
	ScatterPlot chart;
	LineSpec precisionLine, distanceLine, altPrecLine;
	Series precision, distance, altPrec;
	LinearAxis xAxis;
	LinearAxis yAxis;
	DefaultChartHandler handler;
public:
	virtual ~PrecisionRecallChart() {}

	PrecisionRecallChart() :
		labelBox( 0, 0, 0, 0 ),
		precisionLine( FL_GREEN, 2, FL_SOLID ),
		precision( 0, &precisionLine, 0, "precision" ),
		distanceLine( FL_RED, 2, FL_SOLID ),
		distance( 0, &distanceLine, 0, "distance" ),
		altPrecLine( FL_BLUE, 2, FL_DASH ),
		altPrec( 0, &altPrecLine, 0, "alternate precision" ),
		xAxis( 0, 1 ),
		yAxis( 1, 0 ),
		handler(&chart)
		//
	{
		add( &labelBox, Location::West );
		add( &chart, Location::Centre );

		chart.Add( &precision )
			.Add( &distance )
			.Add( &altPrec )
			.SetMargin( 40, 10, 30, 30 )
			.SetAxesCross( 0, 0 )
			.SetXAxis( &xAxis )
			.SetYAxis( &yAxis )
			.SetTickMarks( 10, "%g", 6, "%0.2f" )
			.AddMouseHandler( &handler );

		labelBox.labelfont( FL_COURIER );
		labelBox.labelsize( 14 );
		labelBox.Anchor( 0, 0.5, 5, 0, 0, 0.5 );
	}

	void SetDataSource( uint rankingIdx, const ::Ranking & ranking ) override {
		listItem = &ranking;

		static const Fl_Color paleGoldenrod = fl_rgb_color( 0xee, 0xe8, 0xaa );
		static const Fl_Color cornflowerBlue = fl_rgb_color( 0x64, 0x95, 0xed );
		static const Fl_Color white = fl_rgb_color( 255, 255, 255 );
		static const Fl_Color lightCornflowerBlue = fl_color_average( white, cornflowerBlue, 0.75 );
		static const Fl_Color palette[] = { white, fl_color_average( white, lightCornflowerBlue, 0.5 ), lightCornflowerBlue };
		static const uint paletteSize = sizeof( palette ) / sizeof( palette[0] );

		auto bgColour = palette[rankingIdx % paletteSize];

		chart.SetFillColour( bgColour );

		labelBox.SetLabel( ranking.sequence->Name() );
		labelBox.color( bgColour );

		precision.Clear();
		distance.Clear();
		altPrec.Clear();

		for ( uint i = 0; i < ranking.recall.size(); i++ ) {
			if ( i < ranking.precision.size() ) precision.Add( ranking.recall[i], ranking.precision[i] );
			if ( i < ranking.knn.elements.size() ) distance.Add( ranking.recall[i], ranking.knn.elements[i].first );
		}

		uint numSteps = 100;
		double averagePrecision = 0;
		vector<double> precision;
		SigRank::GetPrecisionRecall(ranking.precision, ranking.recall, numSteps, precision, averagePrecision, rankingIdx );

		for ( uint i = 0; i <= numSteps; i++ ) {
			altPrec.Add( (double) i / numSteps, precision[i] );
		}
	}

	void SetLabelWidth( int labelWidth ) {
		labelBox.size( labelWidth, labelBox.h() );
		size( w(), h() );
	}
};

class PrecisionRecallBrowser : public RankingBrowserBase {
private:
	const int labelWidth = 250;
	Pointers<PrecisionRecallChart> newPointer;

public:

	PrecisionRecallBrowser( int recordHeight = 150 ) : RankingBrowserBase( recordHeight ) {
		SetTitle( "Precision and Distance as function of Recall (ctrl + s to save chart)" );
	}

	virtual ~PrecisionRecallBrowser() {}

	virtual ListItem<::Ranking> * GetItemDisplay() {
		auto p = newPointer();
		p->SetLabelWidth( labelWidth );
		return p;
	}
};
