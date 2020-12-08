#pragma once

#include <FL/Fl.H>
#include <Fl/Fl_Button.H>
#include <Fl/fl_draw.H>
#include <FL/Fl_Int_Input.H>

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
#include <TextDisplay.hpp>
#include "Helpers.hpp"
#include "ScrollingList.hpp"
#include "RankingBrowserBase.hpp"
#include "SigRank.hpp"

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;
using namespace LBGraph;

struct AlignInfo {
	const unordered_map<size_t, SparseSignature *> *sigIndex;
	const vector<Centroid> * vocab;
	uint kmerLength;
	const SimilarityMatrix* matrix;
	Distance threshold;
};

struct Hit {
	double dist;
	const SparseSignature *querySig;
	const SparseSignature *hitSig;
	const AlignInfo *alignInfo;

	Hit(
		double dist,
		const SparseSignature *querySig,
		const SparseSignature *hitSig,
		const AlignInfo *alignInfo
	) :
		dist( dist ),
		querySig( querySig ),
		hitSig( hitSig ),
		alignInfo( alignInfo ) {}
};

class RankingHitView : public virtual ListItem<Hit>, public virtual ScatterPlot::MouseHandler {
	Box numberBox;
	string number;
	Box labelBox;
	string heading;

	Square featureMarker;
	LineSpec seqLine, arcLine;
	Series query, queryFeatures;
	Series hit, hitFeatures;
	Series arcs;
	ScatterPlot chart;
	DefaultChartHandler handler;

public:
	RankingHitView() :
		labelBox( 0, 0, 0, 30 ),
		numberBox( 0, 0, 75, 30 ),
		featureMarker( 0, FL_BLUE, 0 ),
		seqLine( fl_rgb_color( 0x30, 0xff, 0x30 ), 15 ),
		arcLine( fl_rgb_color( 0x30, 0x30, 0xff ), 1 ),
		query( 0, &seqLine, 0 ),
		queryFeatures( &featureMarker, 0, 0 ),
		hit( 0, &seqLine, 0 ),
		hitFeatures( &featureMarker, 0, 0 ),
		arcs( 0, &arcLine, 0 ),
		handler(&chart)
		//
	{
		labelBox.color( FL_WHITE );
		add( &labelBox, BorderLayout::Location::North );

		numberBox.color( FL_WHITE );
		add( &numberBox, BorderLayout::Location::West );

		chart.Add( &query ).Add( &queryFeatures ).Add( &hit ).Add( &hitFeatures ).Add( &arcs )
			.SetMargin( 75, 10, 30, 30 ).SetFillColour( FL_WHITE ).AddMouseHandler(&handler).AddMouseHandler(this);

		add( &chart, BorderLayout::Location::Centre );
	}

	virtual ~RankingHitView() {}

	void Refresh() {
		labelBox.SetLabel( heading );
		labelBox.AutoFit( 3, 3, 3, 3 );
		numberBox.SetLabel( number );
		ResizeImpl();
	}

	virtual void SetDataSource( uint itemIdx, const Hit & item ) override {
		listItem = &item;

		auto &hitSig = *(listItem->hitSig);
		auto &hitSeq = *hitSig.Sequence();
		auto hitLen = hitSeq.Length();

		heading = String::Format( "%f: %s", listItem->dist, hitSeq.DefLine().c_str() );
		number = String::Format( "Hit %d", itemIdx );

		auto &querySig = *(listItem->querySig);
		auto &querySeq = *querySig.Sequence();
		auto queryLen = querySig.Sequence()->Length();

		auto alignInfo = *(listItem->alignInfo);

		query.Clear();
		query.Add( 0, 0.9 );
		query.Add( queryLen, 0.9 );

		vector<SigRank::Feature> queryFeatureList;
		SigRank::GetAllFeatures( querySeq, querySig, alignInfo.kmerLength, *(alignInfo.vocab), alignInfo.threshold, *(alignInfo.matrix), queryFeatureList );

		queryFeatures.Clear();

		for ( auto & f : queryFeatureList ) {
			for ( auto &occurrence : f.occurrence ) {
				queryFeatures.Add( occurrence.kmerPos + alignInfo.kmerLength / 2, 0.9, occurrence.distance );
			}
		}

		hit.Clear();
		hit.Add( 0, 0.1 );
		hit.Add( hitLen, 0.1 );

		vector<SigRank::Feature> hitFeatureList;
		SigRank::GetAllFeatures( hitSeq, hitSig, alignInfo.kmerLength, *(alignInfo.vocab), alignInfo.threshold, *(alignInfo.matrix), hitFeatureList );

		hitFeatures.Clear();

		for ( auto & f : hitFeatureList ) {
			for ( auto &occurrence : f.occurrence ) {
				hitFeatures.Add( occurrence.kmerPos + alignInfo.kmerLength / 2, 0.1, occurrence.distance );
			}
		}

		arcs.Clear();

		uint m = queryFeatureList.size();
		uint n = hitFeatureList.size();
		uint i = 0, j = 0;

		while ( i < m && j < n ) {
			uint x = queryFeatureList[i].centroidPos, y = hitFeatureList[j].centroidPos;

			if ( x < y ) {
				i++;
			}
			else if ( y < x ) {
				j++;
			}
			else {
				for ( auto & qPos : queryFeatureList[i].occurrence ) {
					for ( auto & hPos : hitFeatureList[j].occurrence ) {
						arcs.Add( qPos.kmerPos + alignInfo.kmerLength / 2, 0.9 );
						arcs.Add( hPos.kmerPos + alignInfo.kmerLength / 2, 0.1 );
						arcs.Add( NAN, NAN );
					}
				}
				i++;
				j++;
			}
		}

		chart.XAxis()->SetMax( std::max( queryLen, hitLen ) );
		chart.SetTickMarks( 10, "%0.1f", -1, "" );

		LineSpec line( FL_BLACK, 1 );
		Plus plus( &line, 3 );

		chart.VTics()->Add( 0, 0.9 )->SetLabel( "Query", -10, 0, 1.0, 0.5 ).SetMarker( &plus );
		chart.VTics()->Add( 0, 0.1 )->SetLabel( "Hit", -10, 0, 1.0, 0.5 ).SetMarker( &plus );

		Refresh();
	}

	virtual bool Handle( const LBGraph::ScatterPlot::MouseEvent & event ) override {
		if ( event.source != &chart ) return false;

		if ( event.eventCode == FL_KEYUP ) {
			auto eventKey = Fl::event_key();
			auto eventState = Fl::event_state();

			if ( tolower( eventKey ) == 'd' ) {
				if ( eventState & (FL_CONTROL | FL_COMMAND) ) {
					auto outFile = DefaultChartHandler::ShowDialog();
					SaveAlignmentData( outFile );
					return true;
				}
			}
		}

		return false;
	}

	void SaveAlignmentData( const string & outFile ) {
		ofstream out(outFile);
		CsvWriter w(out);

		auto alignInfo = *(listItem->alignInfo);
		auto &querySig = *(listItem->querySig);
		auto &querySeq = *querySig.Sequence();
		auto queryLen = querySig.Sequence()->Length();
		vector<SigRank::Feature> queryFeatureList;
		SigRank::GetAllFeatures( querySeq, querySig, alignInfo.kmerLength, *(alignInfo.vocab), alignInfo.threshold, *(alignInfo.matrix), queryFeatureList );

		(w << "Feature" << "Kmer Pos" << "Distance" << "Sequence").Ln();
		
		for ( auto & f: queryFeatureList ) {
			for ( auto & occ: f.occurrence) { 
				(w << f.centroidPos << occ.kmerPos << occ.distance << querySeq.IdStr()).Ln();
			}
		}

		auto &hitSig = *(listItem->hitSig);
		auto &hitSeq = *hitSig.Sequence();
		auto hitLen = hitSeq.Length();
		vector<SigRank::Feature> hitFeatureList;
		SigRank::GetAllFeatures( hitSeq, hitSig, alignInfo.kmerLength, *(alignInfo.vocab), alignInfo.threshold, *(alignInfo.matrix), hitFeatureList );
		
		for ( auto & f: hitFeatureList ) {
			for ( auto & occ: f.occurrence) { 
				(w << f.centroidPos << occ.kmerPos << occ.distance << hitSeq.IdStr()).Ln();
			}
		}
	}
};

class RankingDetailListView : public ScrollingList<Hit> {
	Pointers<RankingHitView> newPointer;

public:
	RankingDetailListView( int recordHeight = 200 ) : ScrollingList<Hit>( recordHeight ) {}

	virtual ~RankingDetailListView() {}

	virtual ListItem<Hit> * GetItemDisplay() override {
		return newPointer();
	}

	virtual bool CompareItems( const Hit & lhs, const Hit & rhs ) const override {
		return lhs.dist < rhs.dist;
	}

	virtual string ItemKey( const Hit & listItem ) const override {
		return listItem.hitSig->Sequence()->DefLine();
	}
};

class RankingDetailView : public ListItem<::Ranking> {
private:
	string heading;
	RankingDetailListView content;
	const AlignInfo * alignInfo;
	vector<Hit> hits;

public:
	virtual ~RankingDetailView() {}

	RankingDetailView() {
		content.getListItems = [this]() { return &hits; };
		box( FL_DOWN_FRAME );
		add( &content );
		auto & titleBox = content.TitleBox();
		titleBox.color( FL_WHITE );
		titleBox.labelfont( FL_COURIER );
		titleBox.labelsize( 14 );
	}

	void Init(
		decltype(alignInfo) alignInfo
	) {
		this->alignInfo = alignInfo;
	}

	void SetDataSource( uint idx, const ::Ranking &ranking ) {
		if ( alignInfo == 0 ) throw Exception( "Not correctly initialised!", FileAndLine );

		listItem = &ranking;

		auto queryPos = alignInfo->sigIndex->find( ranking.sequence->Id() );
		auto querySig = queryPos->second;

		hits.clear();
		for ( const auto & element : ranking.knn.elements ) {
			auto hitPos = alignInfo->sigIndex->find( element.second->Id() );
			auto hitSig = hitPos->second;
			hits.emplace_back( element.first, querySig, hitSig, alignInfo );
		}

		heading = "Query Sequence: " + ranking.sequence->DefLine();
		auto & titleBox = content.TitleBox();
		titleBox.SetLabel( heading );
		titleBox.AutoFit( 3, 3, 3, 3 );
		int w, h;
		GetAvailableSize( w, h );
		titleBox.size( w, titleBox.h() );

		content.Clear();
		content.getListItems = [this]() { return &hits; };
		content.GainFocus();
	}

	void Clear() {
		content.Clear();
	}
};

class RankingDetailBrowser : public RankingBrowserBase {
	Pointers<RankingDetailView> newPointer;
	const AlignInfo *alignInfo;

public:
	RankingDetailBrowser() : RankingBrowserBase( 0 ) {
		SetTitle( "Ranking Details Per Query" );
	}

	virtual ~RankingDetailBrowser() {}

	void Init(
		decltype(alignInfo) alignInfo
	) {
		this->alignInfo = alignInfo;
	}

	virtual ListItem<::Ranking> * GetItemDisplay() override {
		if ( alignInfo == 0 ) throw Exception( "Not correctly initialised!", FileAndLine );

		auto p = newPointer();
		p->Init( alignInfo );
		return p;
	}

	virtual void Clear() override {
		for ( int i = 0; i < contentHolder.children(); i++ ) {
			RankingDetailView * item = dynamic_cast<RankingDetailView *>(contentHolder.child( i ));
			item->Clear();
		}

		ScrollingList::Clear();
	}
};
