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

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;
using namespace LBGraph;

class PrecisionDistanceChart : public ListItem<::Ranking> {
private:
	Box labelBox;
	ScatterPlot chart;
	LineSpec precisionLine, recallLine;
	Series precision, recall;
	LinearAxis xAxis;
	LinearAxis yAxis;
	int labelWidth;
public:
	virtual ~PrecisionDistanceChart() {}

	PrecisionDistanceChart() :
		labelWidth( labelWidth ),
		labelBox( 0, 0, labelWidth, 0 ),
		precisionLine( FL_GREEN, 2 ),
		precision( 0, &precisionLine, 0, "precision" ),
		recallLine( FL_BLUE, 2 ),
		recall( 0, &recallLine, 0, "recall" ),
		xAxis( 0, 1 ),
		yAxis( 1, 0 )
		//
	{
		add( &labelBox, Location::West );
		add( &chart, Location::Centre );

		labelBox.labelfont( FL_COURIER );
		labelBox.labelsize( 14 );
		labelBox.anchorX = 0;
		labelBox.anchorY = 0.5;
		labelBox.anchorXto = 0;
		labelBox.anchorYto = 0.5;
		labelBox.offsetX = 5;
		labelBox.offsetY = 0;

		chart.Add( &precision )
			.Add( &recall )
			.SetMargin( 40, 10, 30, 30 )
			.SetAxesCross( 0, 0 )
			.SetXAxis( &xAxis )
			.SetYAxis( &yAxis )
			.SetTickMarks( 10, "%g", 6, "%0.2f" );
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
		recall.Clear();

		for ( uint i = 0; i < ranking.recall.size(); i++ ) {
			if ( i >= ranking.precision.size() ) continue;
			if ( i >= ranking.knn.elements.size() ) continue;

			precision.Add( ranking.knn.elements[i].first, ranking.precision[i] );
			recall.Add( ranking.knn.elements[i].first, ranking.recall[i] );
		}
	}

	void SetLabelWidth( int width ) {
		labelBox.size( width, labelBox.h() );
		size( w(), h() );
	}
};

class PrecisionDistanceBrowser : public RankingBrowserBase {
private:
	const int labelWidth = 250;
	Pointers<PrecisionDistanceChart> newPointer;
public:

	PrecisionDistanceBrowser( int recordHeight = 150 ) : RankingBrowserBase( recordHeight ) {
		SetTitle( "Precision and Recall as function of Distance" );
	}

	virtual ~PrecisionDistanceBrowser() {}

	virtual ListItem<::Ranking> * GetItemDisplay() override {
		auto p = newPointer();
		p->SetLabelWidth( labelWidth );
		return p;
	}
};
