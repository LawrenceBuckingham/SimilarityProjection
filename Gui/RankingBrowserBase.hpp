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
#include "ScrollingList.hpp"

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;
using namespace LBGraph;

class RankingBrowserBase : public ScrollingList<::Ranking> {
public:

	virtual ~RankingBrowserBase() {}

	/**
	 *	Initialise a general-purpose ranking browser, assigning a horizontal row
	 *	to display some kind of results for each query sequence. If rowHeight is
	 *	strictly positive, then the assigned row for each query will be at least
	 *	rowHeight pixels high, with as many records packed in as can fit while
	 *	maintaining this invariant.
	 */
	RankingBrowserBase(int recordHeight = 150 ) : ScrollingList( recordHeight ) {}

	virtual bool CompareItems( const ::Ranking & lhs, const ::Ranking & rhs ) const override {
		auto nameIndex = lhs.sequence->NameIndex();
		auto classIndex = lhs.sequence->ClassIndex();

		if ( nameIndex >= 0 ) {
			auto & lstr = lhs.sequence->Name();
			auto & rstr = rhs.sequence->Name();
			return lstr < rstr;
		}
		else if ( classIndex >= 0 ) {
			auto & lstr = lhs.sequence->Metadata( classIndex );
			auto & rstr = rhs.sequence->Metadata( classIndex );
			return lstr < rstr;
		}
		else {
			return lhs.sequence->Id() < rhs.sequence->Id();
		}
	}

	virtual string ItemKey( const ::Ranking & listItem ) const override {
		return listItem.sequence->DefLine();
	}
};
