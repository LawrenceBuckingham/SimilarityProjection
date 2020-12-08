#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <random>

namespace QutBio {
	struct UniformRealRandom {
		std::mt19937 rng;
		std::uniform_real_distribution<double> distribution;

		UniformRealRandom(unsigned long seed) : rng(seed), distribution(0.0, 1.0) {}

		double operator()() {
			return distribution(rng);
		}

		void Reseed( unsigned long seed ) {
			rng.seed( seed );
		}
	};

	/**
	 *	Uniform integer random values in a sub-range.
	 *
	 *	@tparam	T	Integral numeric type to generate.
	 */
	template<typename T = int>
	struct UniformIntRandom {
		std::mt19937 rng;
		std::uniform_int_distribution<T> distribution;

		/**
		 *	Initialise a UniformIntRandom object.
		 */
		UniformIntRandom(
			/// @param seed	The seed for the underlying random number generator.
			unsigned long seed, 

			/// @param min	The (included) lower bound of the sub-range which 
			///				be sampled by the generator.
			T min = 0,
			
			///	@param max	The (!!! included !!!) upper bound of the sub-range
			///				that will be sampled by the generator.
			T max = 100
		) :
			rng(seed),
			distribution(min, max)
			//
		{}

		/**
		 *	Call without arguments to generate a random value in the _closed_ interval, [min,max].
		 */
		T operator()() {
			auto x = distribution(rng);
			return x;
		}

		/**
		 *	Call with explicit values of min and max to ad-hoc sample a different interval than
		 *	that which is stored.
		 */
		T operator()(T min, T max) {
			std::uniform_int_distribution<T> d(min, max);
			return d(rng);
		}

		/**
		 *	Re-seed the random number generator.
		 */
		void Reseed( 
			/// @param seed	The seed for the underlying random number generator.
			unsigned long seed,

			/// @param min	The (included) lower bound of the sub-range which 
			///				be sampled by the generator.
			T min = 0,

			///	@param max	The (!!! included !!!) upper bound of the sub-range
			///				that will be sampled by the generator.
			T max = 100
		) {
			rng.seed(seed);
			std::uniform_int_distribution<T> newDistribution(min, max);
			distribution = newDistribution;
		}
	};
}
