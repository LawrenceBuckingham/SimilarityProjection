#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <algorithm>
#include <cmath>
#include <functional>
#include <ios>
#include <vector>


#include "Assert.hpp"
#include "CsvIO.hpp"
#include "Distribution.hpp"
#include "Exception.hpp"
#include "NormalDistribution.hpp"

namespace QutBio {
	class ImpossibleDistribution : public Distribution {
	public:
		virtual double Cdf( double t ) { return 0; }
		virtual double Pdf( double t ) { return 0; }
		virtual double InverseCdf( double t ) { return NAN; }
		virtual double Mean( void ) { return NAN; }
		virtual double StdDev( void ) { return NAN; }
	};

	class IntegerDistribution : public Distribution {
	protected:
		int min;
		int max;
		vector<double> p;
		vector<double> f;
		double mu = NAN;
		double sigma = NAN;

	public:
		IntegerDistribution( int min, int max ) : min( min ), max( max ), p( 1 + max - min ), f( 1 + max - min ) {}
		
		IntegerDistribution( int min, int max, const vector<double> & values ) : min( min ), max( max ), p( 1 + max - min ), f( 1 + max - min ) {
			double *p = this->p.data();
			const double *pp = values.data();

			double sum = 0;
			double *f = this->f.data();

			for ( int i = min; i <= max; i++ ) {
				double prob = pp[i - min];

				assert_true( prob >= 0 );

				p[i - min] = prob;
				sum += prob;

				f[i - min] = sum;
			}

			for ( int i = min; i <= max; i++ ) {
				p[i - min] /= sum;
				f[i - min] /= sum;
			}
		}

		template<typename T>
		IntegerDistribution( const Histogram<T> & hist ) {
			min = hist.data.begin()->first;
			max = hist.data.rbegin()->first;
			p.resize( 1 + max - min );
			f.resize( 1 + max - min );

			double *p = this->p.data();

			double sum = 0;
			double *f = this->f.data();

			for ( int i = min; i <= max; i++ ) {
				double prob = hist[i];

				assert_true( prob >= 0 );

				p[i - min] = prob;
				sum += prob;

				f[i - min] = sum;
			}

			for ( int i = min; i <= max; i++ ) {
				p[i - min] /= sum;
				f[i - min] /= sum;
			}
		}

		IntegerDistribution( istream & in ) {
			Parse( in );
		}

		IntegerDistribution( const IntegerDistribution & other ) {
			this->min = other.min;
			this->max = other.max;
			this->p = other.p;
			this->f = other.f;
			this->mu = other.mu;
			this->sigma = other.sigma;
		}

		IntegerDistribution & operator=( const IntegerDistribution & other ) {
			this->min = other.min;
			this->max = other.max;
			this->p = other.p;
			this->f = other.f;
			this->mu = other.mu;
			this->sigma = other.sigma;
			return *this;
		}

		virtual ~IntegerDistribution() {}

		virtual double Cdf( double t ) {
			if ( t < min ) return 0;
			if ( t > max ) return 1;
			return f[(int) floor( t ) - min];
		}

		virtual double Pdf( double t ) {
			if ( t < min ) return 0;
			if ( t > max ) return 0;
			return p[(int) floor( t ) - min];
		}

		virtual double InverseCdf( double t ) {
			if ( t <= 0 ) return min;
			if ( t >= 1 ) return max;

			for ( uint i = 0; i < f.size(); i++ ) {
				if ( f[i] >= t ) {
					return i + min;
				}
			}

			return max;
		}

		virtual double Mean( void ) {
			if ( std::isnan( mu ) ) {
				double *p = this->p.data();
				double sum = 0;

				for ( int i = min; i <= max; i++ ) {
					sum += p[i - min] * i;
				}

				mu = sum;
			}

			return mu;
		}

		virtual double StdDev( void ) {
			if ( std::isnan( sigma ) ) {
				double *p = this->p.data();
				double sumSq = 0;
				double mu = Mean();

				for ( int i = min; i <= max; i++ ) {
					sumSq += p[i - min] * ( i - mu ) * ( i - mu );
				}

				sigma = sqrt( sumSq );
			}

			return sigma;
		}

		double P( int i ) const {
			assert_true( min <= i && i <= max );
			return p[i - min];
		}

		int Min() const { return min; }

		int Max() const { return max; }

		IntegerDistribution Add( const IntegerDistribution & other ) {
			IntegerDistribution &my = *this;
			int newMin = my.min + other.min;
			int newMax = my.max + other.max;

			vector<double> newP( 1 + newMax - newMin );

			for ( int i = min; i <= max; i++ ) {
				for ( int j = other.min; j <= other.max; j++ ) {
					newP[i + j - newMin] += P( i ) * other.P( j );
				}
			}

			return IntegerDistribution( newMin, newMax, newP );
		}

		IntegerDistribution AddParallel( const IntegerDistribution & other ) {
			IntegerDistribution &my = *this;
			int newMin = my.min + other.min;
			int newMax = my.max + other.max;

			vector<double> newP( 1 + newMax - newMin );

#pragma omp parallel
			{
				vector<double> newPPT( 1 + newMax - newMin );

#pragma omp for
				for ( int i = min; i <= max; i++ ) {
					for ( int j = other.min; j <= other.max; j++ ) {
						newPPT[i + j - newMin] += P( i ) * other.P( j );
					}
				}

#pragma omp critical
				{
					for ( int i = 0; i < 1 + newMax - newMin; i++ ) {
						newP[i] += newPPT[i];
					}
				}
			}

			return IntegerDistribution( newMin, newMax, newP );
		}

		/// <summary>
		/// Emits the distribution to standard output in a tabular form.
		/// </summary>
		/// <param name="keyFormat">A function that converts the key to string format.</param>
		/// <param name="keyFormat">A function that converts the key to string format.</param>
		/// <returns>A function that converts the key to string format.</returns>

		ostream & Print(
			ostream & out,
			function<string( double )> valFormat
		) {

			double mu = Mean();
			double sigma = StdDev();
			NormalDistribution norm( mu, sigma );

			double *p = this->p.data();
			out << "x";

			for ( int i = min; i <= max; i++ ) {
				double P = p[i - min];
				double F = f[i - min];

				if ( P > 0 || F > 0 ) {
					out << "\t" << i;
				}
			}

			out << "\nP";

			for ( int i = min; i <= max; i++ ) {
				double P = p[i - min];
				double F = f[i - min];

				if ( P > 0 || F > 0 ) {
					out << "\t" << valFormat( P );
				}
			}
			out << "\nF";

			for ( int i = min; i <= max; i++ ) {
				double P = p[i - min];
				double F = f[i - min];

				if ( P > 0 || F > 0 ) {
					out << "\t" << valFormat( F );
				}
			}

			out << "\nN(" << mu << "," << sigma << ")";

			for ( int i = min; i <= max; i++ ) {
				double P = p[i - min];
				double F = f[i - min];
				double N = norm.Cdf( i );

				if ( P > 0 || F > 0 ) {
					out << "\t" << valFormat( N );
				}
			}

			out << "\nError";

			for ( int i = min; i <= max; i++ ) {
				double P = p[i - min];
				double F = f[i - min];
				double N = norm.Cdf( i );

				if ( P > 0 || F > 0 ) {
					out << "\t" << fabs( N - F );
				}
			}

			out << "\n";

			return out;
		}

		ostream & PrintPdf(
			ostream & out
		) {
			double *p = this->p.data();
			out << "x";

			for ( int i = min; i <= max; i++ ) {
				double P = p[i - min];
				double F = f[i - min];

				if ( P > 0 || F > 0 ) {
					out << "\t" << i;
				}
			}

			out << "\nP";

			for ( int i = min; i <= max; i++ ) {
				double P = p[i - min];
				double F = f[i - min];

				if ( P > 0 || F > 0 ) {
					out << "\t" << setprecision( 17 ) << P;
				}
			}
			out << "\n";
			return out;
		}

		/**
		 *	Consumes two lines for input stream 'in', the first of which
		 *	contains the integer values x_i over which this distribution is
		 *	defined, and the second of which contains the count (or
		 *	probability) of each x_i. The resulting histogram is normalised
		 *	and loaded into the p and f tables.
		 */

	protected:
		void Parse( istream & in ) {
			typedef pair<int, double> P;

			p.clear();
			vector<P> input;

			CsvReader csv( in, '\t' );
			vector<string> r1, r2;

			csv.ReadRecord( r1 );

			if ( r1.size() <= 1 ) {
				r1.clear();
				csv.ReadRecord( r1 );
			}

			csv.ReadRecord( r2 );

			if ( r1.size() != r2.size() || r1.size() <= 1 ) {
				cerr << "r1.size() = " << r1.size() << "\n";
				cerr << "r2.size() = " << r2.size() << "\n";
				cerr << "IntegerDistribution: Unable to parse PDF. Series dimensions do not make sense.\n";
				throw Exception( "IntegerDistribution: Unable to parse PDF. Series dimensions do not make sense.\n", FileAndLine );
			}

			// Field 0 should hold literal 'x' (r1) or 'P' (r2), so skip the first.

			for ( uint i = 1; i < r1.size(); i++ ) {
				int key = atoi( r1[i].c_str() );
				double value = atof( r2[i].c_str() );

				assert_true( value >= 0 );

				input.emplace_back( key, value );
			}

			sort( input.begin(), input.end(), []( P & x, P & y ) { return x.first < y.first; } );

			int newMin = input[0].first;
			int newMax = input.back().first;

			min = newMin;
			max = newMax;
			p.resize( 1 + max - min );
			f.resize( 1 + max - min );

			double sum = 0;

			for ( auto kvp : input ) {
				sum += kvp.second;
				p[kvp.first - min] = kvp.second;
			}

			for ( int i = min; i <= max; i++ ) {
				p[i - min] /= sum;
				f[i - min] = i == 0 ? p[i - min] : f[i-min-1] + p[i - min];
			}
		}

	public:
		/// <summary>
		///	Gets the distribution of the maximum of a subset containing 
		///	subsetSize items from the underlying set.
		/// </summary>

		IntegerDistribution GetMaximum( int subsetSize ) {
			double *f = this->f.data();

			int newMin = max;
			int newMax = min;

			vector<double> Fm( 1 + max - min );

			for ( int i = min; i <= max; i++ ) {
				// http://stats.stackexchange.com/questions/220/how-is-the-minimum-of-a-set-of-random-variables-distributed
				double F = f[i - min];
				Fm[i - min] = pow( F, subsetSize );
			}

			for ( int i = min; i <= max; i++ ) {
				if ( Fm[i - min] > 0 ) {
					newMin = i;
					break;
				}
			}

			for ( int i = max; i >= min; i-- ) {
				if ( Fm[i - min] > 0 ) {
					newMax = i;
					break;
				}
			}

			IntegerDistribution d( newMin, newMax );
			double *fMax = d.f.data();
			double *pMax = d.p.data();

			double fPrev = 0;

			for ( int i = newMin; i <= newMax; i++ ) {
				fMax[i - newMin] = Fm[i - min];
				pMax[i - newMin] = Fm[i - min] - fPrev;
				fPrev = Fm[i - min];
			}

			return d;
		}

		/// <summary>
		///	Gets the distribution of the maximum of a subset containing 
		///	subsetSize items from the underlying set.
		/// </summary>

		IntegerDistribution GetMinimum( int subsetSize ) {
			double *f = this->f.data();

			int newMin = max;
			int newMax = min;

			vector<double> Fm( 1 + max - min );

			for ( int i = min; i <= max; i++ ) {
				// http://stats.stackexchange.com/questions/220/how-is-the-minimum-of-a-set-of-random-variables-distributed
				double F = i == min ? 0 : f[i - 1 - min];
				double Fc = 1 - F;

				if ( F < 1e-10 ) {
					//	Once F is less than about 1e-17, we completely lose the
					//	probability contribution for the left end. To be on the 
					//	safe side I'll stretch the iterative procedure up to 
					//	1e-10.
					double currentTerm = 1;
					int sign = -1;
					Fm[i - min] = 0;

					for ( int j = 1; j <= subsetSize; j++ ) {
						sign = -sign;
						currentTerm *= ( subsetSize - j + 1 ) * F / j;
						Fm[i - min] += sign * currentTerm;
						// cerr << "F = " << F << "\nn = " << subsetSize << "\nj = " << j << "\ncurrentTerm = " << currentTerm << "\nFm = " << Fm[i - min] << "\n";
					}
				}
				else {
					Fm[i - min] = 1 - pow( Fc, subsetSize );
				}
			}

			for ( int i = min; i <= max; i++ ) {
				if ( Fm[i - min] > 0 ) {
					newMin = i;
					break;
				}
			}

			for ( int i = max; i >= min; i-- ) {
				if ( Fm[i - min] > 0 ) {
					newMax = i;
					break;
				}
			}

			IntegerDistribution d( newMin, newMax );
			double *fMax = d.f.data();
			double *pMax = d.p.data();

			double fPrev = 0;

			for ( int i = newMin; i <= newMax; i++ ) {
				fMax[i - newMin] = Fm[i - min];
				pMax[i - newMin] = Fm[i - min] - fPrev;
				fPrev = Fm[i - min];
			}

			return d;
		}

		/// <summary>
		///	Returns a distribution for a conditional event, as determined by the
		///	predicate, which must return true to indicate that the designated value
		///	satisfies the condition.
		/// The arguments passed to the predicate are: x, PDF(x) and CDF(x). 
		///	</summary>

		IntegerDistribution GetConditionalDistribution( function<bool( int, double, double )> predicate ) {
			int newMin = max;
			int newMax = min;

			for ( int i = min; i <= max; i++ ) {
				if ( predicate( i, p[i - min], f[i - min] ) ) {
					newMin = i;
					break;
				}
			}

			for ( int i = max; i >= min; i-- ) {
				if ( predicate( i, p[i - min], f[i - min] ) ) {
					newMax = i;
					break;
				}
			}

			vector<double> newP( 1 + newMax - newMin );

			for ( int i = newMin; i <= newMax; i++ ) {
				if ( predicate( i, p[i - min], f[i - min] ) ) {
					newP[i - newMin] = p[i - min];
				}
			}

			return IntegerDistribution( newMin, newMax, newP );
		}

		void TabulateCdf( vector<double> & x, vector<double> & F ) {
			x.clear();
			F.clear();

			for ( int i = min; i < max; i++ ) {
				double x_i = i;
				double F_i = Cdf( x_i );
				x.push_back( x_i );
				F.push_back( F_i );
			}
		}

		/**
		 *	Computes the theoretical distribution of pairwise k-mer distances
		 *	by propagating pairwise symbol distance distributions through k
		 *
		 *	@tparam <C>	Symbol type, e.g. char, uint8_t
		 *	@tparam <D> Symbol distance type, e.g. Distance
		 *	@returns The distribution of k-mer word distances under the hypothesis 
		 *				that symbols in all kmers are drawn iid from the 
		 *				distribution characterised by symbolHistogram.
		 */
		template<typename C, typename D>
		static IntegerDistribution GetKmerDistanceDistribution(
			/// @param symbolHistogram	Observed frequency distribution of symbols.
			const Histogram<C> &symbolHistogram,

			/// @param symbolDistance Functor which computes pairwise symbol distances.
			function<D( C, C )> & symbolDistance,

			/// @param kmerLength Number of symbols per word: the k in k-mer.
			uint kmerLength
			//
		) {
			Histogram<D> oneMers;
			oneMers.GetOneMerHistogram( symbolHistogram, symbolDistance );

			IntegerDistribution d1( oneMers );
			IntegerDistribution current( d1 );

			for ( uint i = 2; i <= kmerLength; i++ ) {
				IntegerDistribution d = current.Add( d1 );
				current = d;
			}

			return current;
		}

	};
}
