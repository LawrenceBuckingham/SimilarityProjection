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

template<typename T>
class ListItem : public BorderLayout {
public:
	const T * listItem;

	ListItem() : listItem( 0 ) {}

	virtual ~ListItem() {}

	virtual void SetDataSource( uint itemIdx, const T & ranking ) = 0;
};

template<typename T>
class ScrollingList : public BorderLayout, public PropertyChangedEventSource
	//
{
#pragma region General state.
protected:
	const int rowHeight = 25;
	const int buttonWidth = 50;

	GridLayout header;
	Box titleBox;
	FlowLayout buttonHolder;
	VerticalFitLayout contentHolder;
	ActionWidget<Fl_Button> go;
	ActionWidget<Fl_Int_Input> focus;
	Box spacer;
	ValueWidget<string, Fl_Input> searchBox;
	ActionWidget<Fl_Button> searchButton;
	ActionWidget<Fl_Scrollbar> scrollbar;

	const vector<T> * listItems;
	vector<long> index;
	bool sorted = false;

	vector<ListItem<T> *> freeList;
#pragma endregion

public:

	virtual ~ScrollingList() {}

	/**
	 *	Initialise a general-purpose ranking browser, assigning a horizontal row
	 *	to display some kind of results for each query sequence. If rowHeight is
	 *	strictly positive, then the assigned row for each query will be at least
	 *	rowHeight pixels high, with as many records packed in as can fit while
	 *	maintaining this invariant.
	 */
	ScrollingList( int recordHeight = 150 ) :
		PropertyChangedEventSource( this ),
		header( 0, 0, 0, rowHeight * 2 + 5, 2 ),
		titleBox( 0, 0, 0, rowHeight, "Browse Rankings" ),
		buttonHolder( 0, 0, 0, rowHeight ),
		contentHolder( recordHeight ),
		go( 0, 0, buttonWidth, rowHeight, "Go to", [this]() { Refresh(); } ),
		focus( 0, 0, 75, rowHeight, 0, go.Action() ),
		listItems( 0 ),
		spacer( 0, 0, 30, rowHeight ),
		searchBox( 0, 0, 125, rowHeight, "Search ..." ),
		searchButton( 0, 0, 75, rowHeight, "Search", [this]() { Search(); } ),
		scrollbar( 0, 0, 20, 20, 0, [this]() { Scroll(); } )
		//
	{
		add( &header, BorderLayout::Location::North );
		header.add( &titleBox );
		header.add( &buttonHolder );
		add( &contentHolder, BorderLayout::Location::Centre );
		add( &scrollbar, BorderLayout::Location::East );

		scrollbar.type( FL_VERTICAL );

		buttonHolder.SetAlignment( FlowLayout::Align::Centre );
		buttonHolder.add( &focus );
		buttonHolder.add( &go );
		buttonHolder.add( &spacer );
		buttonHolder.add( &searchBox );
		buttonHolder.add( &searchButton );

		SetFocus( 0 );
	}

	void SetItems_defunct( const vector<T> & items ) {
		this->listItems = &items;
		GainFocus();
		Refresh();
	}

	function<decltype(listItems) ()> getListItems;

	void MoveNext() {
		auto n = Focus();
		SetFocus( n + 1 );
	}

	void MovePrev() {
		auto n = Focus();
		SetFocus( n - 1 );
	}

	void MoveFirst() {
		SetFocus( 0 );
	}

	void MoveLast() {
		SetFocus( listItems.size() - contentHolder.Rows() - 1 );
	}

	/**
	 *	Get the index of the current focus of this control (the record being displayed).
	 *	@returns The index of the current focus. This is the value of the focus widget.
	 */
	long Focus() {
		auto v = focus.value();
		long val;

		if ( v && *v ) {
			val = atol( v );
		}
		else {
			val = 0;
			focus.value( "0" );
		}

		auto N = listItems ? (long) listItems->size() : 0;
		auto recordsPerPage = contentHolder.Rows();

		if ( N > recordsPerPage && val > N - recordsPerPage ) {
			val = N - recordsPerPage;
			focus.value( String::Format( "%zu", val ).c_str() );
		}

		return val;
	}

	/**
	 *	Set the index of the record to be displayed in the control.
	 *	@param focus The index of the
	 */
	void SetFocus( long value ) {
		auto N = listItems ? (long) listItems->size() : (long) 0;
		auto recordsPerPage = contentHolder.Rows();

		if ( value > N - recordsPerPage ) value = N - recordsPerPage;
		if ( value < 0 ) value = 0;

		if ( Focus() == value ) return;

		char buff[128];
		snprintf(buff, sizeof(buff), "%ld", value);
		this->focus.value( buff );
		Refresh();
		NotifyPropertyChanged( NAMEOF( Focus ) );
	}

	void Refresh() {
		Sort();

		auto N = listItems ? listItems->size() : 0;

		auto focus = Focus();
		auto recordsPerPage = std::min( (int) N, contentHolder.Rows() );

		int nChildren;

		while ( (nChildren = contentHolder.children()) < recordsPerPage ) {
			ListItem<T> * itemDisplay = GetItemDisplayFreeList();
			contentHolder.add( itemDisplay );
		}

		while ( (nChildren = contentHolder.children()) > recordsPerPage ) {
			ListItem<T> * item = dynamic_cast<ListItem<T> *>(contentHolder.child( nChildren - 1 ));
			freeList.push_back( item );
			contentHolder.remove( nChildren - 1 );
		}

		scrollbar.value( focus );

		if ( recordsPerPage > 0 ) {
			for ( uint i = focus; i < N && i < focus + recordsPerPage; i++ ) {
				ListItem<T> * itemDisplay = dynamic_cast<ListItem<T> *>(contentHolder.child( i - focus ));
				itemDisplay->SetDataSource( i, (*listItems)[index[i]] );
			}
		}

		redraw();
	}

private:
	ListItem<T> * GetItemDisplayFreeList() {
		if ( freeList.empty() ) return GetItemDisplay();
		auto ret = freeList.back();
		freeList.pop_back();
		return ret;
	}

	ListItem<T> * GetItemDisplayFreeList( const T * item ) {
		for ( uint i = 0; i < freeList.size(); i++ ) {
			if ( freeList[i]->listItem == item ) {
				auto ret = freeList[i];
				freeList.erase( freeList.begin() + i );
				return ret;
			}
		}

		return 0;
	}

public:
	virtual ListItem<T> * GetItemDisplay() = 0;
	//{
	//	throw NotImplementedException( FileAndLine );
	//}

	virtual bool CompareItems( const T & lhs, const T & rhs ) const  = 0;
	//{
	//	throw NotImplementedException( FileAndLine );
	//}

	virtual string ItemKey( const T & listItem ) const = 0;
	//{
	//	throw NotImplementedException( FileAndLine );
	//}

private:
	void Sort() {
		if ( sorted || listItems == 0 || listItems->size() == 0 ) return;

		index.resize( listItems->size() );

		for ( long i = 0; i < (long) index.size(); i++ ) {
			index[i] = i;
		}

		std::sort( index.begin(), index.end(), [this]( long lhs, long rhs ) {
			return CompareItems( (*listItems)[lhs], (*listItems)[rhs] );
		} );

		sorted = true;
	}

public:
	void GainFocus() {
		sorted = false;
		focus.value( "0" );
		listItems = getListItems();

		if ( listItems ) {
			scrollbar.bounds( 0, (*listItems).size() + 1 - contentHolder.Rows() );
			scrollbar.slider_size( (*listItems).size() > 0 ? (double) contentHolder.Rows() / (*listItems).size() : 1.0 );
			scrollbar.linesize( 1 );
		}

		Refresh();
	}

private:
	void Search() {
		string key = String::ToLowerCase( String::Trim( searchBox.value() ) );

		if ( key.length() == 0 ) return;

		auto focus = Focus();
		size_t n = listItems->size();

		for ( uint i = focus + 1; i < n; i++ ) {
			string s = String::ToLowerCase( ItemKey( (*listItems)[index[i]] ) );

			if ( s.find( key ) != string::npos ) {
				SetFocus( i );
				return;
			}
		}

		for ( uint i = 0; i < focus; i++ ) {
			string s = String::ToLowerCase( ItemKey( (*listItems)[index[i]] ) );

			if ( s.find( key ) != string::npos ) {
				SetFocus( i );
				return;
			}
		}
	}

public:
	void draw() override {
		BorderLayout::draw();
	}

	void resize( int x, int y, int w, int h ) override {
		BorderLayout::resize( x, y, w, h );
		Refresh();
	}

private:
	void Scroll() {
		SetFocus( scrollbar.value() );
	}

private:
	string title;

public:
	const string & Title() const {
		return title;
	}

	void SetTitle( const string & s ) {
		title = s;
		titleBox.SetLabel( title );
	}

	Box &TitleBox() { return titleBox; }

	/**
	 *	Deletes all cached ListItems and sets items to null.
	 */
	virtual void Clear() {
		listItems = 0;
		freeList.clear();
		Refresh();
	}
};
