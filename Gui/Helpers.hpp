#pragma once

#include <FL/Fl.H>
#include <Index.hpp>
#include <SimpleKmer.hpp>
#include <CsvIO.hpp>
#include <Centroid.hpp>
#include <FastaSequence.hpp>
#include <kNearestNeighbours.hpp>

struct SeqSigIdx { size_t seqIdx, sigIdx; };

struct Ranking {
	const FastaSequence * sequence = 0;
	KnnVector<const FastaSequence *, double> knn;
	vector<double> precision;
	vector<double> recall;

	Ranking(
		const FastaSequence * sequence,
		KnnVector<const FastaSequence *, double> knn,
		vector<double> precision,
		vector<double> recall
	) : sequence( sequence ), knn( knn ), precision( precision ), recall( recall ) {}
};

using pKmerInstance = shared_ptr<const SimpleKmer::Instance>;

class Helpers {
public:
	static void GetCentroid(
		CsvReader & r,
		string & workingStorage,
		size_t kmerLength,
		const LookupTable_<size_t, const FastaSequence> & seqIndex,
		Centroid &c
	) {
		size_t kmerOffset;
		size_t initialClusterSize;
		size_t finalClusterSize;
		size_t finalInstanceCount;
		double purity;
		double entropy;
		r >> workingStorage >> kmerOffset >> initialClusterSize >> finalClusterSize >> finalInstanceCount >> purity >> entropy;

		c.centroid = Helpers::GetKmer( workingStorage, kmerOffset, kmerLength, seqIndex );
		c.initialClusterSize = initialClusterSize;
		c.finalClusterSize = finalClusterSize;
		c.finalInstanceCount = finalInstanceCount;
		c.purity = purity;
		c.entropy = entropy;
	}

	static pKmerInstance GetKmer(
		string & idStr,
		size_t pos,
		size_t kmerLength,
		const LookupTable_<size_t, const FastaSequence> & seqIndex
	) {
		const FastaSequence * seq = 0;

		try {
			auto idNumber = FastaSequence::Register( idStr );
			seq = seqIndex.at( idNumber );
		}
		catch ( out_of_range ex ) {
			throw Exception( "Sequence not found: " + string( idStr ), FileAndLine );
		}

		pKmerInstance ret( new SimpleKmer::Instance( seq, pos ) );
		return ret;
	}

	static void UpdateGui( function<void()> action ) {
		Fl::lock();
		action();
		Fl::unlock();
		Fl::awake();
	}

#define UPDATE_GUI(stmts) { Fl::lock(); {stmts;} Fl::unlock(); Fl::awake();}
};
