#pragma once

#include <FL/Fl.H>
#include <Fl/Fl_Button.H>
#include <FL/Fl_Int_Input.H>

#include <FlowLayout.hpp>
#include <BorderLayout.hpp>
#include <SparseSignature.hpp>
#include <ScrollArea.hpp>

#include <FastaSequence.hpp>
#include <SparseSignature.hpp>
#include <ActionWidget.hpp>

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;

class SignatureCanvas : public Fl_Widget {

#pragma region Focus
private:
	long focus;

public:
	long Focus() {
		return focus;
	}

	void SetFocus( const long val ) {
		if ( val == this->focus ) return;

		this->focus = val;
		redraw();
	}
#pragma endregion

#pragma region Signatures and sequences
private:
	int rowHeight = 25;
	int labelWidth = 0;
	const vector<SparseSignature> * signatures = 0;
	const vector<FastaSequence *> * dbSeqs = 0;
	vector<long> index;
	bool sorted = false;
#pragma endregion

#pragma region Constructor
public:
	SignatureCanvas(
		int rowHeight,
		int labelWidth
	) : Fl_Widget( 0, 0, 0, 0 ), rowHeight( rowHeight ), labelWidth( labelWidth ), focus( 0 ) {}

#pragma endregion

#pragma region RecordPerPage property
	long RecordsPerPage() {
		Fl_Boxtype b = box();
		int hh = h() - Fl::box_dh( b );
		return hh < 0 ? 0 : hh / rowHeight;
	}
#pragma endregion

#pragma region Draw
	void draw() override {
		DrawImpl();
	}

	void DrawImpl() {
		draw_box();
		Fl_Boxtype b = box();
		int xx = x() + Fl::box_dx( b );  // was 9 instead of dx...
		int yy = y() + Fl::box_dy( b );
		int ww = w() - Fl::box_dw( b );
		int hh = h() - Fl::box_dh( b );
		fl_push_clip( xx, yy, ww, hh );

		fl_color( FL_WHITE );
		fl_rectf( xx, yy, ww, hh );

		if ( signatures ) {
			long recordsPerPage = hh / rowHeight;

			int bins = ww - labelWidth - 2;
			vector<double> histogram( bins );

			Sort();

			for ( long i = focus; i < focus + recordsPerPage && i < (long) signatures->size(); i++ ) {
				std::fill( histogram.begin(), histogram.end(), 0 );

				long top = (i - focus) * hh / recordsPerPage;
				long bottom = (i - focus + 1) * hh / recordsPerPage - 1;

				auto & sig = (*signatures)[index[i]];
				auto seq = (*dbSeqs)[index[i]];
				const string & label = nameIndex > 0 ? seq->Metadata( nameIndex ) : seq->IdStr();
				fl_color( FL_BLACK );
				fl_rect( xx + labelWidth, yy + top, ww - labelWidth, bottom - top + 1 );
				fl_draw( label.c_str(), xx, yy + bottom - fl_descent() );

				for ( auto feature : sig ) {
					int from = feature * bins / vocabSize;
					int to = std::max( from, (int) ((feature + 1) * bins / vocabSize - 1) );

					for ( auto j = from; j <= to; j++ ) {
						histogram[j] ++;
					}
				}

				auto max = Util::Max( histogram.begin(), histogram.end(), 0.0 );

				for ( auto &x : histogram ) {
					x = (int) x * 255 / max;
				}

				for ( int i = 0; i < bins; i++ ) {
					fl_color( histogram[i], histogram[i], histogram[i] );
					fl_rect( xx + labelWidth + 1 + i, yy + top + 1, 1, bottom - top - 1 );
				}
			}
		}

		draw_label();
		fl_pop_clip();
	}
#pragma endregion

#pragma region NameIndex Property
private:
	int nameIndex = -1;

public:
	int NameIndex() {
		return nameIndex;
	}

	void SetNameIndex( const int val ) {
		if ( val == this->nameIndex ) return;

		this->nameIndex = val;
		sorted = false;
	}
#pragma endregion

#pragma region ClassIndex Property
private:
	int classIndex = -1;

public:
	int ClassIndex() {
		return classIndex;
	}

	void SetClassIndex( const int val ) {
		if ( val == this->classIndex ) return;

		this->classIndex = val;
		sorted = false;
	}
#pragma endregion

public:
	/// Specify the sequence database which will be displayed. This needs to be 
	///	called after the constructor and before  the first call to draw.
	void SetDbSeqs( decltype(dbSeqs) dbSeqs, decltype(signatures) signatures ) {
		this->dbSeqs = dbSeqs;
		this->signatures = signatures;
		sorted = false;
	}

	void Sort() {
		if ( sorted || dbSeqs == 0 ) return;

		index.resize( dbSeqs->size() );

		for ( long i = 0; i < (long) index.size(); i++ ) {
			index[i] = i;
		}

		if ( nameIndex >= 0 ) {
			std::sort( index.begin(), index.end(), [this]( long lhs, long rhs ) {
				if ( (*dbSeqs)[lhs]->Metadata( nameIndex ) < (*dbSeqs)[rhs]->Metadata( nameIndex ) ) return true;
				return false;
			} );
		}
		else if ( classIndex >= 0 ) {
			std::sort( index.begin(), index.end(), [this]( long lhs, long rhs ) {
				if ( (*dbSeqs)[lhs]->Metadata( classIndex ) < (*dbSeqs)[rhs]->Metadata( classIndex ) ) return true;
				return false;
			} );
		}
		else {
			std::sort( index.begin(), index.end(), [this]( long lhs, long rhs ) {
				if ( (*dbSeqs)[lhs]->Id() < (*dbSeqs)[rhs]->Id() ) return true;
				return false;
			} );
		}

		sorted = true;
	}

private:
	int vocabSize;

public:
	int VocabSize() {
		return vocabSize;
	}

	void SetVocabSize( const int val ) {
		if ( val == this->vocabSize ) return;

		this->vocabSize = val;
	}
};

class SignatureBrowser :
	public BorderLayout
	//
{
#pragma region General state.
	const int rowHeight = 25;
	const int buttonWidth = 50;

	FlowLayout buttonHolder;
	ActionWidget<Fl_Button> first, prev, next, last;
	ActionWidget<Fl_Int_Input> focus;
	Fl_Box label1;

	SignatureCanvas canvas;

	const vector<FastaSequence *> * dbSeqs;
	const vector<SparseSignature> * signatures;
#pragma endregion

public:

	virtual ~SignatureBrowser() {}

#pragma region Constructors
	SignatureBrowser() :
		BorderLayout( 0, 0, 0, 0 ),
		dbSeqs( 0 ),
		signatures( 0 ),
		buttonHolder( 0, 0, 0, rowHeight ),
		label1( 0, 0, 150, rowHeight, "First Sequence:" ),
		first( 0, 0, buttonWidth, rowHeight, "|<", [this]() { MoveFirst(); } ),
		prev( 0, 0, buttonWidth, rowHeight, "<<", [this]() { MovePrev(); } ),
		focus( 0, 0, 75, rowHeight, 0, [this]() { canvas.SetFocus( Focus() ); } ),
		next( 0, 0, buttonWidth, rowHeight, ">>", [this]() { MoveNext(); } ),
		last( 0, 0, buttonWidth, rowHeight, ">|", [this]() { MoveLast(); } ),
		canvas( rowHeight, 250 )
		//
	{
		add( &buttonHolder, BorderLayout::Location::North );
		add( canvas );

		buttonHolder.add( &label1 );
		buttonHolder.add( &first );
		buttonHolder.add( &prev );
		buttonHolder.add( &focus );
		buttonHolder.add( &next );
		buttonHolder.add( &last );

		SetFocus( 0 );

		label1.align( FL_ALIGN_RIGHT | FL_ALIGN_INSIDE );
	}
#pragma endregion 

#pragma region Implementations of button click actions.
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
		SetFocus( signatures->size() - canvas.RecordsPerPage() - 1 );
	}
#pragma endregion

#pragma region Focus property
	/**
	 *	Get the index of the current focus of this control (the record being displayed).
	 *	@returns The index of the current focus. This is the value of the focus widget.
	 */
	long Focus() const {
		auto v = focus.value();
		return v ? (*v ? atol( v ) : -1) : -1;
	}

	/**
	 *	Set the index of the record to be displayed in the control.
	 *	@param focus The index of the
	 */
	void SetFocus( long value ) {
		auto N = signatures ? (long) signatures->size() : 0;
		auto recordsPerPage = canvas.RecordsPerPage();

		if ( value >= N - recordsPerPage ) value = N - recordsPerPage - 1;
		if ( value < 0 ) value = 0;

		if ( Focus() == value ) return;

		this->focus.value( String::Format( "%zu", value ).c_str() );
		canvas.SetFocus( value );
	}
#pragma endregion

#pragma region NameIndex Property
public:
	int NameIndex() {
		return canvas.NameIndex();
	}

	void SetNameIndex( const int val ) {
		canvas.SetNameIndex( val );
	}
#pragma endregion

#pragma region ClassIndex Property
public:
	int ClassIndex() {
		return canvas.ClassIndex();
	}

	void SetClassIndex( const int val ) {
		canvas.SetClassIndex( val );
	}
#pragma endregion

#pragma region DbSeqs Property
public:
	/// Specify the sequence database which will be displayed. This needs to be 
	///	called after the constructor and before  the first call to draw.
	void SetDbSeqs( decltype(dbSeqs) dbSeqs, decltype(signatures) signatures ) {
		this->dbSeqs = dbSeqs;
		this->signatures = signatures;
		canvas.SetDbSeqs( dbSeqs, signatures );
	}
#pragma endregion

#pragma region VocabSize Property
public:
	int VocabSize() {
		return canvas.VocabSize();
	}

	void SetVocabSize( const int val ) {
		canvas.SetVocabSize( val );
	}
#pragma endregion

	void Clear() {
		SetDbSeqs( 0, 0 );
	}
};
