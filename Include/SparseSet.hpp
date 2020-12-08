#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>
#include <string>
#include "Exception.hpp"
#include "Types.hpp"
#include "CsvIO.hpp"

using namespace std;

namespace QutBio {
	class SparseSet :
		public virtual ICsvReader,
		public virtual ICsvWriter {

	protected:
		vector<uint> features;
		bool isOrdered = true;

	public:
		SparseSet() {}

		virtual ~SparseSet() {}

		void Clear() {
			features.clear();
			isOrdered = true;
		}

		void Reserve( uint capacity ) {
			features.reserve( capacity );
		}

		void Add( uint feature ) {
			if ( isOrdered && features.size() > 0 && features.back() == feature ) return;

			isOrdered = isOrdered && ((features.size() == 0) || (features.back() < feature));

			features.push_back( feature );
		}

		size_t Size() const {
			return features.size();
		}

		uint Min() const {
			if ( ! isOrdered ) {
				throw Exception("Set must be ordered!", FileAndLine);
			}

			return features.size() == 0 ? numeric_limits<uint>::max() : features.front();
		}

		uint Max() const {
			if ( ! isOrdered ) {
				throw Exception("Set must be ordered!", FileAndLine);
			}

			return features.size() == 0 ? numeric_limits<uint>::max() : features.back();
		}

		bool Contains( uint feature ) const {
			return (features.size() == 0) ? false
				: (isOrdered) ? BinContains( feature )
				: LinContains( feature );
		}

		vector<uint>::const_iterator begin() const {
			return features.begin();
		}

		vector<uint>::const_iterator  end() const {
			return features.end();
		}

		friend ifstream & operator>>( ifstream &str, SparseSet & set ) {
			set.Clear();
			uint64_t t;
			size_t cardinality = 0;

			str >> cardinality;

			// cerr << "cardinality = " << cardinality << "\n";

			if ( cardinality > 0 ) {
				while ( !(str.eof() || str.peek() == ';') ) {
					str >> t;
					if ( str.fail() )
						break;
					set.Add( t );
					// cerr << "pushed " << t << "\n";
				}

				if ( str.peek() == ';' )
					str.ignore( 1 );
			}
			else {
				char c = str.get();

				while ( c != ';' && !str.eof() ) {
					c = str.get();
				}
			}

			if ( set.features.size() != cardinality ) {
				stringstream s;
				s << "bitSet cardinality " << set.features.size() << " does not match expected value: " << cardinality;
				throw Exception( s.str(), FileAndLine );
			}

			set.Sort();

			return str;
		}

		/**
		 *	Gets the Jaccard similarity between this signature and another.
		 *	Note that this may sort and unique the feature vectors, if they are disordered.
		 *	@param other A SparseSignature against which this signature is to be compared.
		 *	@returns The Jaccard similarity between this signature and the other.
		 */
		double Similarity( const SparseSet & other ) const {
			if ( !isOrdered) {
				throw Exception("This sparse signature must be ordered to allow similarity comparison.", FileAndLine);
			}
			if ( !other.isOrdered ) {
				throw Exception("Other sparse signature must be ordered to allow similarity comparison.", FileAndLine);
			}
			return Jaccard( features, other.features );
		}

		void Write( CsvWriter &w ) const override {
			w << features.size();

			for ( auto idx : features ) w << idx;

			w.Ln();
		}

		void Read( CsvReader &r ) override {
			size_t cardinality;
			r >> cardinality;

			Clear();
			Reserve( cardinality );

			for ( size_t i = 0; i < cardinality && ! r.IsEOL(); i++ ) {
				size_t idx;
				r >> idx;
				Add( idx );
			}

			Sort();
		}

		void Sort() {
			if ( !isOrdered ) {
				std::sort( features.begin(), features.end() );
				std::unique( features.begin(), features.end() );
				isOrdered = true;
			}
		}

		void IntersectWith( const SparseSet & other, SparseSet & intersection ) const {
			if ( ! (isOrdered && other.isOrdered) ) {
				throw Exception( "Sparse sets must be sorted.", FileAndLine );
			}

			intersection.Clear();

			uint m = features.size();
			uint n = other.features.size();
			uint i = 0, j = 0;

			while ( i < m && j < n ) {
				uint x = features[i], y = other.features[j];

				if ( x < y ) {
					i++;
				}
				else if ( y < x ) {
					j++;
				}
				else {
					intersection.Add( x );
					i++;
					j++;
				}
			}
		}

	private:
		bool BinContains( uint feature ) const {
			return std::binary_search( features.begin(), features.end(), feature );
		}

		bool LinContains( uint feature ) const {
			return std::find( features.begin(), features.end(), feature ) != features.end();
		}

		/**
		 *	Gets the Jaccard similarity between two lists of feature numbers.
		 *	@param a A list of feature numbers, arranged in ascending order.
		 *	@param b A list of feature numbers, arranged in ascending order.
		 *	@returns The Jaccard similarity between the two feature lists.
		 */
		static double Jaccard( const vector<uint> & a, const vector<uint> & b ) {
			uint m = a.size();
			uint n = b.size();
			uint i = 0, j = 0, intersect = 0, union_ = 0;

			while ( i < m && j < n ) {
				uint x = a[i], y = b[j];

				union_++;

				if ( x < y ) {
					i++;
				}
				else if ( y < x ) {
					j++;
				}
				else {
					intersect++;
					i++;
					j++;
				}
			}

			union_ += m + n - i - j;

			return (double) intersect / union_;
		}
	};
}
