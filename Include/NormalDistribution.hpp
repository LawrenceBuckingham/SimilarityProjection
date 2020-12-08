#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <cmath>

#include "Distribution.hpp"
#include "Exception.hpp"
#include "Constants.h"
#include "JohnCook.h"

namespace QutBio {
	class NormalDistribution : Distribution {
	private:
		double mu, sigma;
	public:
		NormalDistribution( double mu = 0, double sigma = 1 ) : mu( mu ), sigma( sigma ) {}

		double Cdf( double t ) {
			return Cdf( t, mu, sigma );
		};

		static double Cdf( double t, double mu, double sigma ) {
			return (1 + erf( (t - mu) / (sigma * SQRT_2) )) / 2;
		};

		double Pdf( double t ) {
			// https://en.wikipedia.org/wiki/Normal_distribution
			double twoSigmaSquared = 2 * sigma * sigma;
			double x = t - mu;
			return exp( -x * x / twoSigmaSquared ) / sqrt( M_PI * twoSigmaSquared );
		};

		static double Pdf( double t, double mu, double sigma ) {
			// https://en.wikipedia.org/wiki/Normal_distribution
			double twoSigmaSquared = 2 * sigma * sigma;
			double x = t - mu;
			return exp( -x * x / twoSigmaSquared ) / sqrt( M_PI * twoSigmaSquared );
		};

		double InverseCdf( double p ) {
			double z = JohnCook::JC::NormalCDFInverse( p );
			return z * sigma + mu;
		}

		double Mean( void ) {
			return mu;
		}

		double StdDev( void ) {
			return sigma;
		}

		void SetMean( double mu ) {
			this->mu = mu;
		}

		void SetStdDev( double sigma ) {
			this->sigma = sigma;
		}
	};
}
