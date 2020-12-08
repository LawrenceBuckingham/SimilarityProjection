#pragma once

#include <FL/Fl.H>
#include <Fl/Fl_Button.H>
#include <FL/Fl_Int_Input.H>

#include <FlowLayout.hpp>
#include <BorderLayout.hpp>
#include <TextDisplay.hpp>

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;

class ClusterBrowser :
	public BorderLayout
	//
{
	const int rowHeight = 40;
	FlowLayout buttonHolder;
	Fl_Button first, prev, next, last, up, down, home, end;
	Fl_Int_Input focus, page, recordsPerPage;
	TextDisplay display;
	Fl_Box label1, label2, label3;

	const vector<Centroid> &prototypes;
	const vector<vector<shared_ptr<const SimpleKmer::Instance>>> & clusters;
	function<uint()> getKmerLength;

public:

	virtual ~ClusterBrowser() {}

#pragma region
	ClusterBrowser(
		decltype(prototypes) prototypes,
		decltype(clusters) clusters,
		decltype(getKmerLength) getKmerLength
	) :
		BorderLayout( 0, 0, 0, 0 ),
		buttonHolder( 0, 0, 0, rowHeight ),
		label1( 0, 0, 100, rowHeight, "Cluster:" ),
		first( 0, 0, 40, rowHeight, "|<" ),
		prev( 0, 0, 40, rowHeight, "<<" ),
		focus( 0, 0, 50, rowHeight ),
		next( 0, 0, 40, rowHeight, ">>" ),
		last( 0, 0, 40, rowHeight, ">|" ),
		label2( 0, 0, 100, rowHeight, "K-mer page:" ),
		up( 0, 0, 40, rowHeight, "<<" ),
		down( 0, 0, 40, rowHeight, ">>" ),
		home( 0, 0, 40, rowHeight, "|<" ),
		end( 0, 0, 40, rowHeight, ">|" ),
		page( 0, 0, 50, rowHeight ),
		recordsPerPage( 0, 0, 50, rowHeight ),
		display( 0, 0, 0, 0 ),
		label3( 0, 0, 150, rowHeight, "K-mers Per Page: " ),
		prototypes( prototypes ),
		clusters( clusters ),
		getKmerLength( getKmerLength )
		//
	{
		add( &buttonHolder, BorderLayout::Location::North );
		add( &display, BorderLayout::Location::Centre );

		buttonHolder.add( &label1 );
		buttonHolder.add( &first );
		buttonHolder.add( &prev );
		buttonHolder.add( &focus );
		buttonHolder.add( &next );
		buttonHolder.add( &last );
		buttonHolder.add( &label2 );
		buttonHolder.add( &home );
		buttonHolder.add( &up );
		buttonHolder.add( &page );
		buttonHolder.add( &down );
		buttonHolder.add( &end );
		buttonHolder.add( &label3 );
		buttonHolder.add( &recordsPerPage );
		focus.value( "0" );
		page.value( "1" );
		recordsPerPage.value( "100" );
		display.textfont( FL_COURIER );
		display.textsize( 16 );
		label1.align( FL_ALIGN_RIGHT | FL_ALIGN_INSIDE );
		label2.align( FL_ALIGN_RIGHT | FL_ALIGN_INSIDE );
		label3.align( FL_ALIGN_RIGHT | FL_ALIGN_INSIDE );

#define CALLBACK(widget) widget.callback( widget ## Callback, this)
		CALLBACK( first );
		CALLBACK( next );
		CALLBACK( focus );
		CALLBACK( prev );
		CALLBACK( last );

		CALLBACK( home );
		CALLBACK( end );
		CALLBACK( page );
		CALLBACK( recordsPerPage );
		CALLBACK( down );
		CALLBACK( up );
#undef CALLBACK

		SetFocus( 0 );
	}
#pragma endregion Constructors

#pragma region
	static void nextCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->MoveNext();
	}

	static void prevCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->MovePrev();
	}

	static void firstCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->MoveFirst();
	}

	static void lastCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->MoveLast();
	}

	static void focusCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->Refresh();
	}

	static void homeCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->HomeAction();
	}

	static void endCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->EndAction();
	}

	static void upCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->UpAction();
	}

	static void downCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->DownAction();
	}

	static void pageCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->Refresh();
	}

	static void recordsPerPageCallback( Fl_Widget * sender, void * userData ) {
		auto This = (ClusterBrowser *) userData;
		This->Refresh();
	}

#pragma endregion Static button click event handlers 

#pragma region
	void MoveNext() {
		auto n = Focus();

		if ( n < prototypes.size() - 1 ) {
			SetFocus( n + 1 );
		}
	}

	void MovePrev() {
		auto n = Focus();
		SetFocus( n - 1 );
	}

	void MoveFirst() {
		SetFocus( 0 );
	}

	void MoveLast() {
		SetFocus( prototypes.size() - 1 );
	}

	void HomeAction() {
		SetPage( 1 );
	}

	void EndAction() {
		auto focus = Focus();

		if ( focus > clusters.size() ) return;

		auto & currentKmers = clusters[focus];
		auto recordsPerPage = RecordsPerPage();
		auto numPages = (currentKmers.size() + recordsPerPage - 1) / recordsPerPage;

		SetPage( numPages );
	}

	void UpAction() {
		auto page = Page();

		if ( page > 1 ) SetPage( page - 1 );
	}

	void DownAction() {
		auto page = Page();
		SetPage( page + 1 );
	}
#pragma endregion Implementations of button click actions.

#pragma region Focus property
	/**
	 *	Get the index of the current focus of this control (the record being displayed).
	 *	@returns The index of the current focus. This is the value of the focus widget.
	 */
	size_t Focus() const {
		try {
			return Ulong::Parse( focus.value() );
		}
		catch ( Exception & ex ) {
			return 0;
		}
	}

	/**
	 *	Set the index of the record to be displayed in the control.
	 *	@param focus The index of the
	 */
	void SetFocus( size_t value ) {
		if ( Focus() == value || value >= prototypes.size() ) return;

		this->focus.value( String::Format( "%zu", value ).c_str() );
		Refresh();
	}
#pragma endregion Focus property.

#pragma region Refresh
	void Refresh() {
		uint kmerLength = getKmerLength();

		TRACE;
		display.Clear();

		if ( prototypes.size() == 0 ) return;

		TRACE;
		auto focus = Focus();

		if ( focus > prototypes.size() ) {
			MoveLast();
			return;
		}

		TRACE;
		display( "Cluster %zu:\n", focus );
		auto &proto = prototypes[focus];
		display( "Prototype %s(%zu):"
			"\n%21s: %s"
			"\n%21s: %zu"
			"\n%21s: %zu"
			"\n%21s: %zu"
			"\n%21s: %g"
			"\n%21s: %g"
			"\n\n",
			proto.centroid->sequence->IdStr().c_str(),
			proto.centroid->kmerPosition,
			"Pattern", proto.centroid->Chars( kmerLength ).c_str(),
			"Initial cluster size", proto.initialClusterSize,
			"Final cluster size", proto.finalClusterSize,
			"Final instance count", proto.finalInstanceCount,
			"Class purity", proto.purity,
			"Class entropy", proto.entropy
		);

		TRACE;
		if ( focus > clusters.size() ) return;

		TRACE;
		auto & currentKmers = clusters[focus];
		auto page = Page();
		auto recordsPerPage = RecordsPerPage();
		auto numPages = std::max( (currentKmers.size() + recordsPerPage - 1) / recordsPerPage, (size_t) 1 );

		if ( page > numPages ) {
			page = numPages;
			this->page.value( Ulong::ToString( page ).c_str() );
		}

		TRACE;
		auto begin = (page - 1)*recordsPerPage;
		auto end = page * RecordsPerPage();

		TRACE;
		for ( uint i = begin; i < end; i++ ) {
			if ( i >= currentKmers.size() ) break;

			auto kmer = currentKmers[i];
			display( "%10u\t%s\t%s\n", i, kmer->Chars( kmerLength ).c_str(), kmer->sequence->DefLine().c_str() );
		}
		TRACE;
	}
#pragma endregion Refresh

#pragma region RecordsPerPage
	uint RecordsPerPage() {
		return Uint::Parse( recordsPerPage.value() );
	}

	void SetRecordsPerPage( uint val ) {
		if ( val < 10 ) val = 10;

		auto oldVal = RecordsPerPage();

		if ( oldVal == val ) return;

		recordsPerPage.value( Uint::ToString( val ).c_str() );

		SetPage( 1 );
	}
#pragma endregion RecordsPerPage

#pragma region Page property
	uint Page() {
		return Uint::Parse( page.value() );
	}

	void SetPage( uint val ) {
		TRACE;
		if ( prototypes.size() == 0 ) return;

		TRACE;
		auto focus = Focus();
		auto &current = clusters[focus];

		if ( current.size() == 0 ) val = 1;

		auto recordsPerPage = RecordsPerPage();
		auto numPages = std::max( (size_t) 1, (current.size() + recordsPerPage - 1) / recordsPerPage );

		if ( val > numPages ) val = numPages;

		TRACE;
		if ( Page() == val ) return;

		TRACE;
		page.value( Int::ToString( val ).c_str() );
		TRACE;
		Refresh();
	}
#pragma endregion Page property
};
