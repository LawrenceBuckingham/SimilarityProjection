#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <thread>
#include <SimilarityMatrix.hpp>
#include <SimpleKmer.hpp>
#include <Centroid.hpp>
#include <SparseSignature.hpp>
#include <FastaSequence.hpp>
#include <kNearestNeighbours.hpp>
#include "Helpers.hpp"
#include <Exception.hpp>

class SigRank {
public:

	/**
	 *	Initialise a ranking set ready for RankSignatures.
	 *
	 *	@param[in] database A list of sparse binary signatures, some (or possibly all) of
	 *		which will be selected as query points.
	 *
	 *	@param[in] queryIndices
	 *
	 *	@param[in] maxMatches The number of results required for each query.
	 *
	 *	@param[out] rankings A list of ::Ranking objects which will be prepared
	 *		for RankSignatures.
	 *
	 *	@post rankings.size = queryIndices.size and
	 *		for all i in dom rankings,
	 *			rankings[i].sequence = database[i].sequence and
	 *			rankings[i].knn.capacity = maxMatches and
	 *			rankings[i].precision.size = maxMatches and
	 *			rankings[i].recall.size = maxMatches
	 */
	static void Setup(
		const vector<SparseSignature> & database,
		const vector<size_t> & queryIndices,
		int maxMatches,
		vector<::Ranking> & rankings
	) {
		rankings.clear();

		KnnVector<const FastaSequence *, double> knn( maxMatches, -1 );
		vector<double> precision;
		vector<double> recall;

		for ( auto i : queryIndices ) {
			auto & sig = database[i];
			rankings.emplace_back( sig.Sequence(), knn, precision, recall );
			auto & ranking = rankings.back();

			auto &precision = ranking.precision;
			precision.clear();
			precision.reserve( maxMatches );

			auto &recall = ranking.recall;
			recall.clear();
			recall.reserve( maxMatches );
		}
	}

	/**
	 *	Identify the k-nearest signatures in the training set (the members of
	 *	which are indexed in a posting list) and rank them in ascending order
	 *	of distance from the designated query. This is intended to be called
	 *	from the body of a parallel loop, so all data structures are required
	 *	to be initialised with necessary capacity before entry.
	 *
	 *	The function will also compute interpolated precision and recall values
	 *	at each recovered document if class labels are available.
	 *
	 *	@param[in] querySig The query signature.
	 *
	 *	@param[in] signatures A list of sparse binary feature vectors.
	 *		Let (s == signatures[i]) be a signature.
	 *		Then (j in s.elements) iff (sequence i contains feature j).
	 *
	 *	@param[in] sigPostingList An inverted index mapping feature numbers to the
	 *		sequences that contain them.
	 *		Let f be a feature.
	 *		Then (i in sigPostingList[f]) iff (f in signatures[i].elements).
	 *
	 *	@param[in] classPostingList An inverted index mapping class labels to the sequences
	 *		that are annotated with those labels.
	 *		Let c be a class.
	 *		Then (i in classPostingList[c]) iff (c in signatures[i].sequence->classes).
	 *
	 *	@param[in/out] processed A bit-set which has been previously initialised
	 *		with capacity to record all values in [0,signatures.size()). Used for
	 *		working storage.
	 *
	 *	@param[in/out] ranking A ::Ranking object with pre-allocated members sufficient
	 *		to store the required number of ranked results, together with any precision
	 *		and recall values that may be computed. The number of results returned
	 *		is equal to the capacity of ranking.knn.
	 */
	static void RankSignatures(
		const SparseSignature & querySig,
		const vector<SparseSignature> & signatures,
		const vector<vector<size_t>> & sigPostingList,
		const unordered_map<size_t, vector<size_t>> & classPostingList,
		BitSet & processed,
		::Ranking & ranking
	) {
		processed.Clear();
		ranking.knn.clear();

		for ( auto c : querySig ) {
			for ( size_t d : sigPostingList[c] ) {
				if ( !processed.Contains( d ) ) {
					processed.Insert( d );
					auto &dbSig = signatures[d];
					double distance = 1.0 - querySig.Similarity( dbSig );

					if ( ranking.knn.canPush( distance ) ) {
						ranking.knn.push( dbSig.Sequence(), distance );
					}
				}
			}
		}

		ranking.knn.sort();
	}

	static uint CountRelevant(
		const FastaSequence * querySeq,
		const unordered_map<size_t, vector<size_t>> & classPostingList,
		BitSet & processed
	) {
		processed.Clear();
		uint result = 0;

		for ( auto classId : querySeq->Classes() ) {
			auto iter = classPostingList.find( classId );

			if ( iter == classPostingList.end() ) {
				// Some classes simply will not exist in the training set, e.g. singleton representative of genus/species.
				continue;
			}

			auto & relatedIndices = iter->second;

			for ( auto i : relatedIndices ) {
				if ( !processed.Contains( i ) ) {
					processed.Insert( i );
					result++;
				}
			}
		}

		return result;
	}

	/**
	 *	Compute interpolated precision and recall values at each recovered document
	 *	if class labels are available.
	 *
	 *	@param[in] classPostingList An inverted index mapping class labels to the sequences
	 *		that are annotated with those labels.
	 *		Let c be a class.
	 *		Then (i in classPostingList[c]) iff (c in signatures[i].sequence->classes).
	 *
	 *	@param[in/out] ranking A ::Ranking object with fully populated k-nearest
	 */
	static void ComputePrecisionRecall(
		const unordered_map<size_t, vector<size_t>> & classPostingList,
		BitSet & processed,
		::Ranking & ranking
	) {
		auto querySeq = ranking.sequence;
		auto &precision = ranking.precision;
		auto &recall = ranking.recall;

		precision.clear();
		recall.clear();

		if ( querySeq->ClassIndex() >= 0 ) {
			uint relevant = 0,
				numRetrieved = ranking.knn.elements.size(),
				totalRelevant = CountRelevant( querySeq, classPostingList, processed );

			for ( uint i = 0; i < numRetrieved; i++ ) {
				auto hitSeq = ranking.knn.elements[i].second;

				if ( querySeq->IsRelated( hitSeq ) ) {
					relevant++;
				}

				precision.push_back( (double) relevant / (i + 1) );
				recall.push_back( totalRelevant == 0 ? 1 : (double) relevant / totalRelevant );
			}

			for ( uint i = 1; i < numRetrieved; i++ ) {
				uint j = numRetrieved - i - 1;

				if ( precision[j] < precision[j + 1] ) {
					precision[j] = precision[j + 1];
				}
			}
		}
	}

	static void ComputePrecisionRecall(
		size_t dbSize,
		const unordered_map<size_t, vector<size_t>> & classPostingList,
		vector<::Ranking> & rankings
	) {
		BitSet processed( dbSize );

		for ( auto & ranking : rankings ) {
			auto numRetrieved = ranking.knn.elements.size();

			if ( ranking.precision.size() != numRetrieved || ranking.recall.size() != numRetrieved ) {
				SigRank::ComputePrecisionRecall( classPostingList, processed, ranking );
			}
		}
	}

	static void GenerateSignature(
		const FastaSequence &seq,
		const size_t kmerLength,
		const vector<Centroid> & vocab,
		const Distance threshold,
		const SimilarityMatrix & matrix,
		SparseSignature & sig
	) {
		sig.SetSequence( &seq );

		const auto x = sig.Sequence()->Sequence().data();
		const auto l = sig.Sequence()->KmerCount( kmerLength );
		const auto M = vocab.size();

		for ( size_t m = 0; m < M; m++ ) {
			auto y = vocab[m].centroid->Bytes();

			for ( size_t j = 0; j < l; j++ ) {
				const auto dist = matrix.Difference( x + j, y, kmerLength );

				if ( dist <= threshold ) {
					sig.Add( m );
					break;
				}
			}
		}
	}

	static void GenerateSignatureDigrams(
		const FastaSequence &seq,
		const size_t kmerLength,
		const vector<Centroid> & vocab,
		const Distance threshold,
		const SimilarityMatrix & matrix,
		SparseSignature & sig
	) {
		sig.Clear();
		sig.SetSequence( &seq );

		auto & digrams = sig.Sequence()->Digrams();
		const auto x = digrams.data();
		const auto l = seq.KmerCount(kmerLength);
		const auto M = vocab.size();

		for ( size_t m = 0; m < M; m++ ) {
			auto y = vocab[m].centroid->Digrams();

			for ( size_t j = 0; j < l; j++ ) {
				const auto dist = matrix.DigramDifference( x + j, y, kmerLength / 2 );

				if ( dist <= threshold ) {
					sig.Add( m );
					break;
				}
			}
		}
	}

	/// Information about location of a feature on a sequence.
	struct Feature {
		struct Occurrence {
			/// Zero origin index of the matching k-mer in sequence.
			size_t kmerPos;

			/// Distance from k-mer to centroid.
			Distance distance;

			Occurrence( size_t kmerPos = 0, Distance distance = 0 ) : kmerPos( kmerPos ), distance( distance ) {}
		};

		/// Zero-origin position of centroid in vocabulary.
		size_t centroidPos;

		/// Zero origin index of the matching k-mer in sequence.
		vector<Occurrence> occurrence;

		Feature( size_t centroidPos = 0 ) : centroidPos( centroidPos ) {}
	};

	/**
	 *	Gets a list containing the annotated features present in
	 */
	static void GetAllFeatures(
		const FastaSequence & seq,
		const SparseSet & sig,
		const size_t kmerLength,
		const vector<Centroid> & vocab,
		const Distance threshold,
		const SimilarityMatrix & matrix,
		vector<SigRank::Feature> & features
	) {
		features.clear();
		const auto x = seq.Sequence().data();
		const auto l = seq.KmerCount( kmerLength );

		for ( auto centroidPos : sig ) {
			auto y = vocab[centroidPos].centroid->Bytes();

			for ( size_t kmerPos = 0; kmerPos < l; kmerPos++ ) {
				const auto dist = matrix.Difference( x + kmerPos, y, kmerLength );

				if ( dist <= threshold ) {
					if ( features.size() == 0 || features.back().centroidPos != centroidPos ) {
						features.emplace_back( centroidPos );
					}

					features.back().occurrence.emplace_back( kmerPos, dist );
				}
			}
		}
	}

	static void GetPrecisionRecall(
		const vector<double> & precision,
		const vector<double> & recall,
		const uint numSteps,
		vector<double> & prec,
		double & averagePrecision,
		size_t recordNumber
	) {
		if ( precision.size() != recall.size() ) {
			auto message = String::Format(
				"precision and recall lists must be non-empty and have the same size.\n"
				"record number = %zu\n"
				"precision.size() = %zu\n"
				"   recall.size() = %zu",
				recordNumber, precision.size(), recall.size()
			);
			throw Exception( message, FileAndLine );
		}

		prec.clear();

		double sum = 0;
		uint n = 0;
		uint i = 0;

		auto r = [numSteps]( uint i ) { return double( i ) / numSteps; };

		if ( precision.size() > 0 ) {
			double lastPrecision = precision[0];

			for ( uint j = 0; j < recall.size(); j++ ) {
				if ( precision[j] == lastPrecision ) {
					while ( i <= numSteps && r( i ) <= recall[j] ) {
						prec.push_back( lastPrecision );
						i++;
					}
				}
				else {
					lastPrecision = precision[j];
				}
			}
		}

		while ( i <= numSteps ) {
			prec.push_back( 0 );
			i++;
		}

		for ( auto p : prec ) {
			sum += p;
			n++;
		}

		averagePrecision = sum / n;
	}

	static void GetPrecisionRecall(
		const vector<::Ranking> & rankings,
		const uint numSteps,
		vector<vector<double>> & prec,
		vector<double> & averagePrecision,
		double & meanAveragePrecision
	) {
		prec.resize( rankings.size() );
		averagePrecision.resize( rankings.size() );
		double sum = 0;
		int n = 0;

		for ( uint i = 0; i < rankings.size(); i++ ) {
			GetPrecisionRecall( rankings[i].precision, rankings[i].recall, numSteps, prec[i], averagePrecision[i], i );
			sum += averagePrecision[i];
		}

		meanAveragePrecision = rankings.size() > 0 ? sum / rankings.size() : 0;
	}

	static void Save( vector<::Ranking> & rankings, string &fileName ) {
		ofstream f( fileName );
		CsvWriter w( f );

		uint n = rankings.size();

		w << "rankings" << n << '\n';
		for ( auto & ranking : rankings ) {
			SaveRanking( w, ranking );
		}

		if ( std::any_of( rankings.begin(), rankings.end(), []( const ::Ranking & r ) { return r.precision.size() > 0 || r.recall.size() > 0; } ) ) {
			w << '\n' << "precision" << '\n';

			for ( auto & ranking : rankings ) {
				SaveVector( w, ranking.sequence->IdStr(), ranking.precision );
			}

			w << '\n' << "recall" << '\n';

			for ( auto & ranking : rankings ) {
				SaveVector( w, ranking.sequence->IdStr(), ranking.recall );
			}
		}
	}

	static void SaveVector( QutBio::CsvWriter &w, const string & idStr, const vector<double> & vec ) {
		w << idStr << vec.size();

		for ( auto x : vec ) {
			w << x;
		}

		w << '\n';
	}

	static void SaveRanking( QutBio::CsvWriter &w, ::Ranking & ranking ) {
		w << ranking.sequence->IdStr() << ranking.knn.capacity;

		for ( auto & neighbour : ranking.knn ) {
			w << neighbour.second->IdStr() << neighbour.first;
		}

		w.Ln();
	}

	static void Load(
		const string &fileName,
		const LookupTable_<size_t, const FastaSequence> & dbIndex,
		const unordered_map<size_t, vector<size_t>> * classPostingList,
		function<void( const string &tag )> errorMsg,
		function<void( size_t maxMatches )> setMaxMatches,
		vector<::Ranking> & rankings
	) {
		ifstream f( fileName );
		CsvReader r( f );
		string tag, queryId, hitId;
		size_t n, m = 0, maxMatches = 0;
		double d;

		r >> tag >> n;

		if ( tag != "rankings" ) {
			errorMsg( tag );
			return;
		}

		rankings.clear();

		unordered_map<size_t, size_t> index;

		for ( size_t i = 0; i < n && !r.IsEOF(); i++ ) {
			r >> queryId >> m;

			if ( m > maxMatches ) maxMatches = m;

			auto seqId = FastaSequence::Register( queryId );
			auto iter = dbIndex.find( seqId );

			if ( iter == dbIndex.end() ) {
				cerr << "Query key not found in signature index '" + queryId << "'. No more records will be parsed.\n";
				break;
			}

			auto querySig = iter->second;

			index[seqId] = rankings.size();
			KnnVector<const FastaSequence *, double> knn_( m, -1 );
			vector<double> precision;
			vector<double> recall;
			rankings.emplace_back( querySig, knn_, precision, recall );
			auto & currentRanking = rankings.back();
			KnnVector<const FastaSequence *, double> &knn = currentRanking.knn;

			while ( !r.IsEOL() ) {
				r >> hitId >> d;
				auto seqId = FastaSequence::Register( hitId );
				auto iter = dbIndex.find( seqId );

				if ( iter == dbIndex.end() ) {
					throw Exception( "Hit key not found in signature index: " + queryId, FileAndLine );
				}

				auto hitSig = iter->second;
				knn.push( hitSig, d );
			}
		}

		setMaxMatches( maxMatches );

		tag.resize( 0 );

		while ( !r.IsEOF() ) {
			r >> tag;

			if ( tag == "precision" ) break;
		}

		if ( tag == "precision" ) {
			for ( size_t i = 0; i < n && !r.IsEOF(); i++ ) {
				r >> queryId >> m;

				auto seqId = FastaSequence::Register( queryId );
				auto iter = index.find( seqId );

				if ( iter == index.end() ) {
					throw Exception( "Query key not found in triple index: " + queryId, FileAndLine );
				}

				::Ranking & ranking = rankings[(*iter).second];
				vector<double> &prec = ranking.precision;
				prec.clear();
				prec.resize( m, 0 );

				for ( uint i = 0; i < m && !r.IsEOL(); i++ ) {
					double p;
					r >> p;
					prec[i] = p;
				}
			}

			while ( !r.IsEOF() ) {
				r >> tag;

				if ( tag == "recall" ) break;
			}

			if ( tag == "recall" ) {
				for ( size_t i = 0; i < n && !r.IsEOF(); i++ ) {
					r >> queryId >> m;

					auto seqId = FastaSequence::Register( queryId );
					auto iter = index.find( seqId );

					if ( iter == index.end() ) {
						throw Exception( "Query key not found in triple index: " + queryId, FileAndLine );
					}

					::Ranking & ranking = rankings[(*iter).second];
					vector<double> &recall = ranking.recall;
					recall.clear();
					recall.resize( m, 0 );

					for ( uint i = 0; i < m && !r.IsEOL(); i++ ) {
						double p;
						r >> p;
						recall[i] = p;
					}
				}
			}
		}
		else if (classPostingList) {
			SigRank::ComputePrecisionRecall( dbIndex.size(), *classPostingList, rankings );
		}

	}

};

struct PrecisionRecallStats {
	uint numSteps = 100;
	vector<vector<double>> prec;
	vector<double> averagePrecision;
	double meanAveragePrecision = 0;

	PrecisionRecallStats & Clear() {
		numSteps = 100;
		prec.resize( 0 );
		averagePrecision.resize( 0 );
		meanAveragePrecision = 0;
		return *this;
	}

	PrecisionRecallStats & Update( const vector<::Ranking> & rankings, uint numSteps ) {
		const uint N = rankings.size();

		if ( N == 0 ) {
			Clear();
			return *this;
		}

		prec.resize( N );
		averagePrecision.resize( N );
		meanAveragePrecision = 0;
		SigRank::GetPrecisionRecall( rankings, numSteps, prec, averagePrecision, meanAveragePrecision );
		return *this;
	}
};
