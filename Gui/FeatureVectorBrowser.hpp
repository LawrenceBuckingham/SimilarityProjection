#pragma once

#include <FL/Fl.H>
#include <Fl/Fl_Button.H>
#include <FL/Fl_Int_Input.H>

#include <FlowLayout.hpp>
#include <BorderLayout.hpp>
#include <SparseSignature.hpp>
#include <ScrollArea.hpp>

#include <FastaSequence.hpp>
#include <SparseFeatureVector.hpp>
#include <ActionWidget.hpp>

#include "ScrollingList.hpp"
#include <Pointers.hpp>

#undef TRON
// #define TRON 1
#include <db.hpp>

using namespace LBFL;

class FeatureVectorBrowser : public ScrollingList<SparseFeatureVector>
	//
{
	struct Canvas : public Fl_Widget {
		size_t * vocabSize;
		const SparseFeatureVector * featureVector = 0;
		vector<double> histogram;

		Canvas (size_t * vocabSize): Fl_Widget(10, 10, 10, 10, 0), vocabSize(vocabSize) {}

		virtual ~Canvas() {}

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

			int bins = ww - 2;
			histogram.resize( bins );

			std::fill( histogram.begin(), histogram.end(), 0 );

			long top = 0;
			long bottom = hh - 1;

			auto & sig = *featureVector;
			fl_color( FL_BLACK );
			fl_rect( xx, yy + top, ww, bottom - top + 1 );

			for ( auto feature : sig ) {
				int from = feature.key * bins / (*vocabSize);
				int to = std::max( from, (int) ((feature.key + 1) * bins / (*vocabSize) - 1) );

				for ( auto j = from; j <= to; j++ ) {
					histogram[j] += feature.weight;
				}
			}

			auto max = Util::Max( histogram.begin(), histogram.end(), 0.0 );

			for ( auto &x : histogram ) {
				x = (int) x * 255 / max;
			}

			for ( int i = 0; i < bins; i++ ) {
				fl_color( histogram[i], histogram[i], histogram[i] );
				fl_rect( xx + 1 + i, yy + top + 1, 1, bottom - top - 1 );
			}

			fl_pop_clip();
		}
	};

	struct Item : public ListItem<SparseFeatureVector> {
		Box labelBox;
		Canvas canvas;

		Item(size_t * vocabSize): labelBox( 0, 0, 250, 0 ), canvas(vocabSize) {
			labelBox.SetFillColour(FL_WHITE).Anchor( 0, 0.5, 3, 0, 0, 0.5 );
			add( &labelBox, LBFL::BorderLayout::Location::West );
			add( &canvas,  LBFL::BorderLayout::Location::Centre );
		}

		virtual void SetDataSource( uint itemIdx, const SparseFeatureVector & featureVector )  {
			listItem = &featureVector;
			canvas.featureVector = &featureVector;
			labelBox.SetLabel(featureVector.Sequence()->Name());
		}

	};

	const vector<SparseFeatureVector> & featureVectors;
	size_t & vocabSize;

	const vector<FastaSequence *> * dbSeqs;

public:
	virtual ~FeatureVectorBrowser() {}

	FeatureVectorBrowser(
		decltype(featureVectors) featureVectors,
		size_t & vocabSize,
		int recordHeight = 25
		//
	) :
		ScrollingList(recordHeight),
		dbSeqs( dbSeqs ),
		featureVectors( featureVectors ),
		vocabSize(vocabSize)
		//
	{
		getListItems = [this]() { return &(this->featureVectors); };
	}

	void SetDbSeqs( decltype(dbSeqs) value ) {
		dbSeqs = value;
	}

	Pointers<Item> newItem;

	virtual ListItem<SparseFeatureVector> * GetItemDisplay() {
		auto item = newItem(&vocabSize);
		return item;
	}

	virtual bool CompareItems( const SparseFeatureVector & lhs, const SparseFeatureVector & rhs ) const  {
		auto lSeq = lhs.Sequence();
		auto rSeq = rhs.Sequence();
		auto nameIndex = lSeq->NameIndex();
		auto classIndex = lSeq->ClassIndex();

		if ( nameIndex >= 0 ) {
			auto & lstr = lSeq->Name();
			auto & rstr = rSeq->Name();
			return lstr < rstr;
		}
		else if ( classIndex >= 0 ) {
			auto & lstr = lSeq->Metadata( classIndex );
			auto & rstr = rSeq->Metadata( classIndex );
			return lstr < rstr;
		}
		else {
			return lSeq->Id() < rSeq->Id();
		}
	}

	virtual string ItemKey( const SparseFeatureVector & listItem ) const {
		return  listItem.Sequence()->DefLine();
	}

};
