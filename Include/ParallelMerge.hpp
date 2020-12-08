#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cstdlib>
#include <thread>

namespace QutBio {
	class ParallelMerge {
	public:
		template<typename T>
		static void sort(T * src, size_t n, T *dst, int(*cmp)(const void * x, const void * y), size_t qsort_limit) {
#pragma omp parallel
			{
#pragma omp single
				sort2(src, 0, n, dst, cmp, qsort_limit);
			}
		}

	private:
		template<typename T>
		static void merge(T * X, size_t start, size_t n, T * tmp, int(*cmp)(const void * x, const void * y)) {
			(cerr << "merge " << start << ".." << (start + n - 1) << endl).flush();

			size_t i = 0;
			size_t j = n / 2;
			size_t ti = 0;

			while ( i < n / 2 && j < n ) {
				if ( cmp(X[start + i], X[start + j]) < 0 ) {
					tmp[start + ti] = X[start + i];
					ti++; i++;
				}
				else {
					tmp[start + ti] = X[start + j];
					ti++; j++;
				}
			}
			while ( i < n / 2 ) { /* finish up lower half */
				tmp[start + ti] = X[start + i];
				ti++; i++;
			}
			while ( j < n ) { /* finish up upper half */
				tmp[start + ti] = X[start + j];
				ti++; j++;
			}
			memcpy(X + start, tmp + start, n*sizeof(T));

		} // end of merge()

		template<typename T>
		static void sort2(T * X, size_t start, size_t n, T * tmp, int(*cmp)(const void * x, const void * y), size_t qsort_limit) {
			if ( n <= qsort_limit ) {
				qsort(X, n, sizeof(T), cmp);
				return;
			}

#pragma omp sections
			{
#pragma omp section
				sort2(X, start, (n / 2), tmp, cmp, qsort_limit);

#pragma omp section
				sort2(X, start + (n / 2), n - (n / 2), tmp, cmp, qsort_limit);
			}

			merge(X, start, n, tmp, cmp);
		}

	};
}

