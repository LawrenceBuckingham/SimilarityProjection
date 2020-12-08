#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <string>
#include <mutex>
#include "Array.hpp"
#include "EnumBase.hpp"

using namespace QutBio;

namespace QutBio {
	class FragmentAggregationMode : public EnumBase {
	private:
		FragmentAggregationMode( string literal, int value ) : EnumBase( literal, value ) {}
		static mutex m;

	public:
		static FragmentAggregationMode * BestOfBest() {
			unique_lock < mutex > lck{m};
			static FragmentAggregationMode value("BestOfBest", 0);
			return &value;
		}

		static FragmentAggregationMode * Hausdorff() {
			unique_lock < mutex > lck{m};
			static FragmentAggregationMode value("Hausdorff", 1);
			return &value;
		}

		static FragmentAggregationMode * HausdorffAverage() {
			unique_lock < mutex > lck{ m };
			static FragmentAggregationMode value( "HausdorffAverage", 2 );
			return &value;
		}

		static FragmentAggregationMode * HausdorffAverageAverage() {
			unique_lock < mutex > lck{ m };
			static FragmentAggregationMode value( "HausdorffAverageAverage", 3 );
			return &value;
		}

		static FragmentAggregationMode * Slice() {
			unique_lock < mutex > lck{m};
			static FragmentAggregationMode value("Slice", 4);
			return &value;
		}

		static FragmentAggregationMode * SliceVertical() {
			unique_lock < mutex > lck{m};
			static FragmentAggregationMode value("SliceVertical", 5);
			return &value;
		}

		static FragmentAggregationMode * SliceNoFollow() {
			unique_lock < mutex > lck{m};
			static FragmentAggregationMode value("SliceNoFollow", 6);
			return &value;
		}

		static FragmentAggregationMode * SliceVerticalNoFollow() {
			unique_lock < mutex > lck{m};
			static FragmentAggregationMode value("SliceVerticalNoFollow", 7);
			return &value;
		}

		static FragmentAggregationMode * Projector() {
			unique_lock < mutex > lck{m};
#pragma warning( disable : 4640 )
			static FragmentAggregationMode value("Projector", 8);
#pragma warning( default : 4640 )
			return &value;
		}

		static FragmentAggregationMode * ProjectorAMP() {
			unique_lock < mutex > lck{m};
#pragma warning( disable : 4640 )
			static FragmentAggregationMode value("ProjectorAMP", 9);
#pragma warning( default : 4640 )
			return &value;
		}

		static FragmentAggregationMode * ProjectorBitEmbedding() {
			unique_lock < mutex > lck{m};
#pragma warning( disable : 4640 )
			static FragmentAggregationMode value("ProjectorBitEmbedding", 10);
#pragma warning( default : 4640 )
			return &value;
		}

		static FragmentAggregationMode * ProjectorSlice() {
			unique_lock < mutex > lck{m};
#pragma warning( disable : 4640 )
			static FragmentAggregationMode value("ProjectorSlice", 11);
#pragma warning( default : 4640 )
			return &value;
		}

		static vector<EnumBase*> Values () {
			vector<EnumBase*> values {
				BestOfBest(), 
				Hausdorff(), 
				HausdorffAverage(),
				HausdorffAverageAverage(),
				Slice(),
				SliceVertical(),
				SliceNoFollow(),
				SliceVerticalNoFollow(),
				Projector(),
				ProjectorAMP(),
				ProjectorBitEmbedding(),
				ProjectorSlice(),
				};
			return values;
		}
	};
}
