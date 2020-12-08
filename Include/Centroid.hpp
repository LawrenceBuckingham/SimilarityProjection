#pragma once

#include "SimpleKmer.hpp"
#include "CsvIO.hpp"

#define UNLIKELY (-1)

struct Centroid :
	public virtual ICsvWriter {
	shared_ptr<const SimpleKmer::Instance> centroid = 0;
	size_t initialClusterSize = 0;
	size_t finalClusterSize = 0;
	size_t finalInstanceCount = 0;
	double purity = UNLIKELY;
	double entropy = UNLIKELY;

	Centroid() {}

	Centroid(
		shared_ptr<const SimpleKmer::Instance> centroid,
		size_t initialClusterSize = 0,
		size_t finalClusterSize = 0,
		size_t finalInstanceCount = 0,
		double purity = UNLIKELY,
		double entropy = UNLIKELY
	) :
		centroid( centroid ),
		initialClusterSize( initialClusterSize ),
		finalClusterSize( finalClusterSize ),
		finalInstanceCount( finalInstanceCount ),
		purity( purity ),
		entropy( entropy ) {}

	virtual ~Centroid() {}

	void Write( CsvWriter & w ) const override {
		w << (*centroid)
			<< initialClusterSize
			<< finalClusterSize
			<< finalInstanceCount
			<< purity
			<< entropy;
	}

	static bool LessByInitSize( const Centroid & lhs, const Centroid & rhs ) {
		return lhs.initialClusterSize < rhs.initialClusterSize;
	}

	static bool LessByFinalSize( const Centroid & lhs, const Centroid & rhs ) {
		return lhs.finalClusterSize < rhs.finalClusterSize;
	}

	static bool GreaterByInitSize( const Centroid & lhs, const Centroid & rhs ) {
		return lhs.initialClusterSize > rhs.initialClusterSize;
	}

	static bool GreaterByFinalSize( const Centroid & lhs, const Centroid & rhs ) {
		return lhs.finalClusterSize == rhs.finalClusterSize ? lhs.initialClusterSize > rhs.initialClusterSize : lhs.finalClusterSize > rhs.finalClusterSize;
	}

	friend bool operator==( const Centroid & lhs, const Centroid & rhs ) {
		return lhs.centroid == rhs.centroid
			&& lhs.entropy == rhs.entropy
			&& lhs.finalClusterSize == rhs.finalClusterSize
			&& lhs.finalInstanceCount == rhs.finalInstanceCount
			&& lhs.initialClusterSize == rhs.initialClusterSize
			;
	}

	friend bool operator!=( const Centroid & lhs, const Centroid & rhs ) {
		return !operator==( lhs, rhs );
	}
};
