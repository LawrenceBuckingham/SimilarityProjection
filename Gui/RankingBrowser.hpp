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
#include <Pointers.hpp>

#include <FastaSequence.hpp>
#include <SparseSignature.hpp>
#include <ActionWidget.hpp>
#include <Box.hpp>
#include <ScatterPlot.hpp>
#include <VerticalFitLayout.hpp>
#include "Helpers.hpp"
#include "RankingBrowserBase.hpp"
#include "DefaultChartHandler.hpp"

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;
using namespace LBGraph;

class RankingChart : public ListItem<::Ranking> {
protected:
	Box labelBox;
	ScatterPlot chart;
	LineSpec precisionLine, recallLine, distanceLine;
	Series precision, recall, distance;
	LinearAxis xAxis;
	LinearAxis yAxis;
	int labelWidth;
	DefaultChartHandler handler;
public:
	virtual ~RankingChart() {}

	RankingChart( int labelWidth = 200 ) :
		labelWidth( labelWidth ),
		labelBox( 0, 0, labelWidth, 0 ),
		precisionLine( FL_GREEN, 2, FL_SOLID ),
		recallLine( FL_BLUE, 2, FL_SOLID ),
		distanceLine( FL_RED, 2, FL_SOLID ),
		precision( 0, &precisionLine, 0, "precision" ),
		recall( 0, &recallLine, 0, "recall" ),
		distance( 0, &distanceLine, 0, "distance" ),
		xAxis( 0, 1 ),
		yAxis( 1, 0 ),
		handler(&chart)
		//
	{
		add( &labelBox, Location::West );
		add( &chart, Location::Centre );

		labelBox.labelfont( FL_COURIER );
		labelBox.labelsize( 14 );
		labelBox.Anchor( 0, 0.5, 5, 0, 0, 0.5 );

		chart.Add( &precision )
			.Add( &recall )
			.Add( &distance )
			.SetMargin( 40, 10, 30, 30 )
			.SetAxesCross( 0, 0 )
			.SetXAxis( &xAxis )
			.SetYAxis( &yAxis )
			.AddMouseHandler(&handler);
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

		labelBox.SetLabel( ranking.sequence->Name() );
		labelBox.color( bgColour );

		chart.SetFillColour( bgColour );

		precision.Clear();

		for ( uint i = 0; i < ranking.precision.size(); i++ ) {
			double y = ranking.precision[i];
			precision.Add( i + 1, y );
		}

		recall.Clear();

		for ( uint i = 0; i < ranking.recall.size(); i++ ) {
			double y = ranking.recall[i];
			recall.Add( i + 1, y );
		}

		distance.Clear();

		for ( uint i = 0; i < ranking.knn.elements.size(); i++ ) {
			double y = ranking.knn.elements[i].first;
			distance.Add( i + 1, y );
		}

		xAxis.SetMax( ranking.precision.size() );
		chart.SetTickMarks( 10, "%g", 6, "%0.2f" );
	};

	void SetLabelWidth( int width ) {
		labelBox.size( width, labelBox.h() );
		size( w(), h() );
	}
};

class RankingBrowser : public RankingBrowserBase {
	const int labelWidth = 250;
	Pointers< RankingChart> newPointer;

public:

	RankingBrowser( int recordHeight = 150 ) : RankingBrowserBase( recordHeight ) {
		SetTitle( "Precision, Recall, and Distance as function of Rank" );
	}

	virtual ~RankingBrowser() {}

	virtual ListItem<::Ranking> * GetItemDisplay() {
		auto p = newPointer();
		p->SetLabelWidth( labelWidth );
		return p;
	}
};
