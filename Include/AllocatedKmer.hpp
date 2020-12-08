#pragma once
#include <atomic>

#include "SimilarityMatrix.hpp"
#include "Kmer.hpp"
#include "KmerCluster.hpp"

namespace QutBio {
#if USE_OMP
#include <omp.h>
	class AllocatedKmer : public Kmer {
		bool isAllocated{ 0 };
		omp_lock_t lock;
	public:
		AllocatedKmer(const Substring & substring) : Kmer(substring) {
			omp_init_lock(&lock);
		}

		AllocatedKmer(const AllocatedKmer & other) : Kmer(other) {
			omp_init_lock(&lock);
			isAllocated = other.isAllocated;
		}

		void operator=(const AllocatedKmer & other) {
			Kmer::operator=(other);
			omp_init_lock(&lock);
			isAllocated = other.isAllocated;
		}

		virtual ~AllocatedKmer() {
			omp_destroy_lock(&lock);
		}

		bool IsAllocated() {
			return isAllocated;
		}

		void Allocate() {
			isAllocated = true;
		}

		void Lock() {
			omp_set_lock(&lock);
		}

		void Unlock() {
			omp_unset_lock(&lock);
		}
	};
#else
	struct AllocatedKmer : Kmer {
		void * connection;

		AllocatedKmer(const Substring & substring) : Kmer(substring) {}

		virtual ~AllocatedKmer() {}

		void Lock() {}

		void Unlock() {}
	};
#endif

}
