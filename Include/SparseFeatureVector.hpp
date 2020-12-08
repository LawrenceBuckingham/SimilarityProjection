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
	/**
	 *	A sparse feature bag data structure; suitable for bag-of-words applications.
	 */
	class SparseFeatureVector :
		public virtual ICsvReader,
		public virtual ICsvWriter
		//
	{
	public:
		struct Feature {
			size_t key;
			double weight;

			Feature( decltype(key) key, decltype(weight) value ) : key( key ), weight( value ) {}

			friend bool operator<( const Feature & lhs, const Feature & rhs ) {
				return lhs.key < rhs.key;
			}

			friend bool operator==( const Feature & lhs, const Feature & rhs ) {
				return lhs.key == rhs.key;
			}

			friend bool operator!=( const Feature & lhs, const Feature & rhs ) {
				return lhs.key == rhs.key;
			}

			struct Hash {
				size_t operator()( const Feature & val ) const noexcept {
					return _Hash_impl::hash( val.key );
				}
			};
		};

	protected:
		const FastaSequence * sequence;
		vector<Feature> features;
		bool isOrdered = true;

	public:
		SparseFeatureVector( const FastaSequence * sequence = 0 ) : sequence( sequence ) {}

		virtual ~SparseFeatureVector() {}

		void Clear() {
			features.clear();
			isOrdered = false;
		}

		void Reserve( uint capacity ) {
			features.reserve( capacity );
		}

		void Add( size_t key, double weight = 1 ) {
			if ( isOrdered && features.size() > 0 && features.back().key == key ) features.back().weight += weight;

			if ( isOrdered ) {
				if ( features.size() == 0 || features.back().key < key ) {
					features.emplace_back( key, weight );
				}
				else {
					auto feature = BinFind( key );

					if ( feature == features.end() || feature->key > key ) {
						features.emplace_back( key, weight );
						isOrdered = false;
					}
					else {
						double &featureWeight = feature->weight;
						featureWeight += weight;
					}
				}
			}
			else {
				auto feature = LinFind( key );

				if ( feature == features.end() ) {
					features.emplace_back( key, weight );
				}
				else {
					double &featureWeight = feature->weight;
					featureWeight += weight;
				}
			}
		}

		decltype(sequence) Sequence() const {
			return sequence;
		}

		void SetSequence( decltype(sequence) value ) {
			sequence = value;
		}

		size_t Size() const {
			return features.size();
		}

		bool Contains( decltype(Feature::key) key ) const {
			return (features.size() == 0) ? false
				: (isOrdered) ? BinFind( key ) != features.end()
				: LinFind( key ) != features.end();
		}

		decltype(Feature::key) MaxKey() const {
			return (features.size() == 0) ? 0
				: (isOrdered) ? features.back().key
				: Util::Max( features.begin(), features.end(), []( const Feature & feature ) { return feature.key; }, 0 );
		}

		decltype(Feature::key) MinKey() const {
			return (features.size() == 0) ? 0
				: (isOrdered) ? features.front().key
				: Util::Min( features.begin(), features.end(), []( const Feature & feature ) { return feature.key; }, 0 );
		}

		/**
		 *	Get a reference to the weight of a designated feature.
		 *	@param key The feature number to dereference.
		 *	@throws out_of_range If the feature is not present, out_of_range is thrown.
		 */
		double & Weight( size_t key ) const {
			if ( features.size() == 0 ) throw out_of_range( "key not found" );
			auto feature = isOrdered ? BinFind( key ) : LinFind( key );
			if ( feature == features.end() ) throw out_of_range( "key not found" );
			return feature->weight;
		}

		auto begin() const {
			return features.begin();
		}

		auto end() const {
			return features.end();
		}

		/**
		 *	Gets the scalar product between this signature and another.
		 *	Note that this may sort and unique the feature vectors, if they are disordered.
		 *	@param other A SparseFeatureVector against which this signature is to be compared.
		 *	@returns The scalar product between this signature and the other.
		 */
		double Dot( const SparseFeatureVector & other ) const {
			if ( !isOrdered ) {
				throw Exception( "This sparse feature vector must be ordered to allow operations with other objects.", FileAndLine );
			}
			if ( !other.isOrdered ) {
				throw Exception( "Other sparse feature vector must be ordered to allow operations with other objects.", FileAndLine );
			}

			return Dot( features, other.features );
		}

		void Write( CsvWriter &w ) const override {
			w << features.size();

			for ( auto idx : features ) w << idx.key << idx.weight;

			w.Ln();
		}

		void Read( CsvReader &r ) override {
			size_t cardinality;
			r >> cardinality;

			Clear();
			Reserve( cardinality );

			for ( size_t i = 0; i < cardinality && !r.IsEOL(); i++ ) {
				size_t key;
				double weight;
				r >> key >> weight;
				Add( key, weight );
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

	private:
		static bool FeatureLessThanKey( const Feature & feature, size_t key ) {
			return feature.key < key;
		}

		decltype(features)::iterator BinFind( decltype(Feature::key) feature ) const {
			auto iter = std::lower_bound( features.begin(), features.end(), feature, FeatureLessThanKey );

			// Yuck! Convert const iterator to iterator, so I can use the return value to update the count.
			vector<Feature> & aaarg = (vector<Feature> &)features;
			return aaarg.erase( iter, iter );
		}

		decltype(features)::iterator  LinFind( decltype(Feature::key) feature ) const {
			auto begin = features.begin();
			auto end = features.end();

			for ( auto i = begin; i != end; i++ ) {
				if ( i->key == feature ) {
					// Yuck! Convert const iterator to iterator, so I can use the return value to update the count.
					vector<Feature> & aaarg = (vector<Feature> &)features;
					return aaarg.erase( i, i );
				}
			}

			// Yuck! Convert const iterator to iterator, so I can use the return value to update the count.
			vector<Feature> & aaarg = (vector<Feature> &)features;
			return aaarg.erase( end, end );
		}

		/**
		 *	Gets the Scalar product between two lists of feature numbers.
		 *	@param a A list of features, arranged in ascending order by key.
		 *	@param b A list of features, arranged in ascending order by key.
		 *	@returns The Jaccard similarity between the two feature lists.
		 */
		static double Dot( const vector<Feature> & a, const vector<Feature> & b ) {
			size_t m = a.size();
			size_t n = b.size();
			size_t i = 0, j = 0;
			double product = 0;

			while ( i < m && j < n ) {
				auto x = a[i].key;
				auto y = b[j].key;

				if ( x < y ) {
					i++;
				}
				else if ( y < x ) {
					j++;
				}
				else {
					product += a[i].weight * b[j].weight;
					i++;
					j++;
				}
			}

			return product;
		}

		/**
		 *	Gets the Scalar product between two lists of feature numbers.
		 *	@param a A list of features, arranged in ascending order by key.
		 *	@param b A list of features, arranged in ascending order by key.
		 *	@returns The Jaccard similarity between the two feature lists.
		 */
		static double Distance( const vector<Feature> & a, const vector<Feature> & b ) {
			size_t m = a.size();
			size_t n = b.size();
			size_t i = 0, j = 0;
			double distance = 0;

			while ( i < m && j < n ) {
				auto x = a[i].key;
				auto y = b[j].key;

				double t;

				if ( x < y ) {
					t = b[j].weight;
					i++;
				}
				else if ( y < x ) {
					t = a[i].weight;
					j++;
				}
				else {
					t = a[i].weight - b[j].weight;
					i++;
					j++;
				}

				distance += t * t;
			}

			return sqrt( distance );
		}

		/**
		 *	Gets the Jaccard similarity between two lists of feature numbers.
		 *	@param a A list of feature numbers, arranged in ascending order.
		 *	@param b A list of feature numbers, arranged in ascending order.
		 *	@returns The Jaccard similarity between the two feature lists.
		 */
		static double Jaccard( const vector<Feature> & a, const vector<Feature> & b ) {
			uint m = a.size();
			uint n = b.size();
			uint i = 0, j = 0, intersect = 0, union_ = 0;

			while ( i < m && j < n ) {
				uint x = a[i].key, y = b[j].key;

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

	public:

		/**
		 *	Derive a posting list from sparse signature.
		 *	@pre (dbSigs.size() == M)
		 *		 and ((i in dom(selectedSignatures)) => (0 <= selectedSignatures[i] < M))
		 *	@post ((i in dom(index) and (j in ran(index[i]))) => (i in ran(dbSigs[j].indices)))
		 */

		static void CreatePostingList(
			const vector<SparseFeatureVector> &dbSigs,
			const vector<size_t> &selectedSignatures,
			vector<vector<size_t>> & index
		) {
			const uint D = selectedSignatures.size();

			uint max_index = 0;

			for ( auto d : selectedSignatures ) {
				for ( const auto & i : dbSigs[d].features ) {
					if ( i.key > max_index ) max_index = i.key;
				}
			}

			index.clear();
			index.resize( max_index + 1 );

			for ( auto & list: index ) {
				list.resize(0);
			}

			for ( auto d : selectedSignatures ) {
				for ( auto i : dbSigs[d].features ) {
					index[i.key].push_back( d );
				}
			}
		}

		/**
		 *	Derive a posting list from sparse signature.
		 *	@pre (dbSigs.size() == M)
		 *		 and ((i in dom(selectedSignatures)) => (0 <= selectedSignatures[i] < M))
		 *	@post ((i in dom(index) and (j in ran(index[i]))) => (i in ran(dbSigs[j].indices)))
		 */

		static void CreatePostingList(
			const vector<SparseFeatureVector> &dbSigs,
			const vector<size_t> &selectedSignatures,
			unordered_map<size_t, vector<size_t>> & index
		) {
			const uint D = selectedSignatures.size();
			index.clear();

			for ( auto d : selectedSignatures ) {
				for ( auto i : dbSigs[d].features ) {
					index[i.key].push_back( d );
				}
			}
		}

	};
}
