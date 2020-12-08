#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

namespace QutBio {
	class Distribution {
	public:
		virtual ~Distribution() {}

		virtual double Cdf( double t ) = 0;
		virtual double Pdf( double t ) = 0;
		virtual double InverseCdf( double t ) = 0;
		virtual double Mean( void ) = 0;
		virtual double StdDev( void ) = 0;

		// Gets (min, max) such that min = inf{t:Cdf(t) > 0} and
		// max = sup{t:Cdf(t) < 1}.
		void GetSupport( double & min, double & max, double epsilon = 1e-10 ) {
			double hi, lo;

			// Get min:
			hi = Mean();
			lo = hi - 100 * StdDev();

			while ( fabs(hi - lo) > epsilon ) {
				min = ( lo + hi ) / 2;
				//( cerr << "lo = " << lo << "\n" ).flush();
				//( cerr << "hi = " << hi << "\n" ).flush();
				//( cerr << "fabs(hi - lo) = " << fabs( hi - lo ) << "\n" ).flush();
				//( cerr << "min = " << min << "\n" ).flush();

				if ( Cdf( min ) <= 0 ) {
					lo = min;
				}
				else {
					hi = min;
				}
			}

			min = ( lo + hi ) / 2;

			// Get max:
			lo = Mean();
			hi = lo + 100 * StdDev();

			while ( fabs( hi - lo ) > epsilon ) {
				max = ( lo + hi ) / 2;
				//(cerr << "max = " << max << "\n").flush();

				if ( Cdf( max ) >= 1 ) {
					hi = max;
				}
				else {
					lo = max;
				}
			}

			max = ( lo + hi ) / 2;
		}
	};

	class ScaledDistribution : public Distribution {
		double scale;
		Distribution &baseDistribution;
	public:
		ScaledDistribution(
			double scale,
			Distribution &baseDistribution
			//
		) : scale( scale ), baseDistribution( baseDistribution ) {}

		virtual ~ScaledDistribution() {}

		virtual double Cdf( double t ) {
			return baseDistribution.Cdf( t / scale );
		}

		virtual double Pdf( double t ) {
			return baseDistribution.Pdf( t / scale );
		}

		virtual double InverseCdf( double t ) {
			return baseDistribution.InverseCdf(t) * scale;
		}

		virtual double Mean( void ) {
			return scale * baseDistribution.Mean();
		}

		virtual double StdDev( void ) {
			return fabs(scale) * baseDistribution.StdDev();
		}

	};

}
