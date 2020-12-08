#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <map>
#include <vector>
#include <iostream>

#include "Distribution.hpp"
#include "Function.hpp"
#include "Histogram.hpp"

namespace QutBio {

	class DiscreteDistribution : public Distribution {
	private:
		Histogram<double> pmf;
		Histogram<double> cdf;
		double mu;
		double sigma;
	public:
		Histogram<double> const & Pmf() { return pmf; }
		Histogram<double> const & Cdf() { return cdf; }

		/// <summary>
		/// Emits the distribution to standard output in a tabular form.
		/// </summary>
		void Print(
			ostream & out,
			Func1<double, string> keyFormat,
			Func1<double, string> valFormat
			) {
			out << "x_i" << "," << "p_i" << "," << "F_i" << endl;

			for ( auto k : pmf.data ) {
				double key = k.first;
				double p = k.second;
				double f = cdf.data[key];
				out << keyFormat(key) << "," << valFormat(p) << "," << valFormat(f) << endl;
			}

			out << endl << "mu," << valFormat(mu) << endl << "sigma," << valFormat(sigma) << endl;
		}

		void SetCdf(Histogram<double> & cdf_) {
			cdf.data.clear();
			pmf.data.clear();

			double maxCdf = 0;
			bool isFirst = true;

			for ( auto k : cdf_.data ) {
				double key = k.first;
				double value = k.second;

				cdf.data[key] = value;

				if ( isFirst ) {
					pmf.data[key] = value;
					isFirst = false;
				}
				else {
					pmf.data[key] = value - maxCdf;
				}

				if ( value > maxCdf ) {
					maxCdf = value;
				}
			}

			double sum_p_x = 0;
			double sum_p_x_2 = 0;

			for ( auto k : cdf_.data ) {
				double key = k.first;
				cdf.data[key] /= maxCdf;
				pmf.data[key] /= maxCdf;
				double p_x = key * pmf.data[key];
				double p_x_2 = p_x * key;
				sum_p_x += p_x;
				sum_p_x_2 += p_x_2;
				//Console.Error( string.Format( "{0}, {1}, {2}, {3}", key, this.pmf[key], p_x, p_x_2 ) );
			}

			mu = sum_p_x;

			//Console.Error( string.Format( "mu = {0}, mu^2 = {1}, sum_p_x_2 = {2}", this.mu, mu * mu, sum_p_x_2 ) );

			sigma = sqrt(sum_p_x_2 - mu * mu);
		}

		void SetPmf(Histogram<double> cdf_) {
			cdf.data.clear();
			pmf.data.clear();

			double maxCdf = 0;
			bool isFirst = true;

			for ( auto k : cdf_.data ) {
				double key = k.first;
				double value = k.second;

				pmf.data[key] = value;

				maxCdf += value;

				if ( isFirst ) {
					cdf.data[key] = value;
					isFirst = false;
				}
				else {
					cdf.data[key] = maxCdf;
				}
			}

			double sum_p_x = 0;
			double sum_p_x_2 = 0;

			for ( auto k : cdf_.data ) {
				double key = k.first;
				cdf.data[key] /= maxCdf;
				pmf.data[key] /= maxCdf;
				double p_x = key * pmf.data[key];
				double p_x_2 = p_x * key;
				sum_p_x += p_x;
				sum_p_x_2 += p_x_2;
				//Console.Error( string.Format( "{0}, {1}, {2}, {3}", key, this.pmf[key], p_x, p_x_2 ) );
			}

			mu = sum_p_x;

			//Console.Error( string.Format( "mu = {0}, mu^2 = {1}, sum_p_x_2 = {2}", this.mu, mu * mu, sum_p_x_2 ) );

			sigma = sqrt(sum_p_x_2 - mu * mu);
		}

		/// <summary>
		///	Gets the distribution of the minimum of a subset containing 
		///	subsetSize items from the underlying set.
		/// </summary>

		void GetMinimumDistribution(
			int subsetSize,
			DiscreteDistribution & outDist
			) {
			Histogram<double> max_cdf;

			for ( auto k : cdf.data ) {
				// http://stats.stackexchange.com/questions/220/how-is-the-minimum-of-a-set-of-random-variables-distributed
				double x = k.first;
				double F = k.second;
				double Fc = 1 - F;
				double Fm = 1 - pow(Fc, subsetSize);
				max_cdf.data[x] = Fm;

				if ( Fm == 1 ) break;
			}

			outDist.SetCdf(max_cdf);
		}

		/// <summary>
		///	Gets the distribution of the maximum of a subset containing 
		///	subsetSize items from the underlying set.
		/// </summary>

		void GetMaximumDistribution(
			int subsetSize,
			DiscreteDistribution & outDist
			) {
			Histogram<double> max_cdf;

			for ( auto k : cdf.data ) {
				// http://stats.stackexchange.com/questions/220/how-is-the-minimum-of-a-set-of-random-variables-distributed
				double x = k.first;
				double F = k.second;
				double Fm = pow(F, subsetSize);
				max_cdf.data[x] = Fm;

				if ( Fm == 1 ) break;
			}

			outDist.SetCdf(max_cdf);
		}

		double Cdf(double t) {
			if ( t < cdf.data.begin()->first ) {
				return 0;
			}
			else if ( t >= cdf.data.rbegin()->first ) {
				return 1;
			}

			//  Reference: lower_bound returns iterator to first element with key 
			//  greater than or equal to t.
			auto low = cdf.data.lower_bound(t);

			if ( low->first == t ) {
				return low->second;
			}
			else if ( low == cdf.data.end() ) {
				return 1;
			}
			else if ( low == cdf.data.begin() ) {
				return 0;
			}
			else {
				auto prev = low;
				--prev;
				double x0 = prev->first, y0 = prev->second;
				double x1 = low->first, y1 = low->second;
				assert_true(x0 < t && x1 > t);
				return y0 + (t - x0) * (y1 - y0) / (x1 - x0);
			}
		}

		double Pdf(double t) {
			// TODO: Review this arbitrary quantity.
			const double delta = 1e-5;
			return (Cdf(t + delta) - Cdf(t - delta)) / (2 * delta);
		}

		double InverseCdf(double t) {
			// TODO: consider a better solution based on directly traversing 
			//  the tree, if that is possible. 
			if ( t < 0 || t > 1 ) {
				return NAN;
			}

			vector<pair<double,double>> vec;
			
			for ( auto & rec : cdf.data ){
				vec.push_back( rec );
			}

			const double epsilon = 1e-10;
			double lo = vec.begin()->first;
			double hi = vec.rbegin()->first;
			double mid = (lo + hi) / 2;

			while ( hi - lo >= epsilon ) {
				mid = (lo + hi) / 2;
				double f = Cdf(mid);
				if ( f > t ) {
					hi = mid;
				}
				else if ( f < t ) {
					lo = mid;
				}
				else {
					break;
				}
			}

			return mid;
		}

		double Mean(void) {
			return mu;
		}

		double StdDev(void) {
			return sigma;
		}

		/// <summary>
		///     Eliminates unnecessary entries for which p == 0.
		/// </summary>

		void Cleanup(void) {
			vector<double> toRemove;

			for ( auto kvp : pmf.data ) {
				if ( kvp.second <= 0 ) {
					toRemove.push_back(kvp.first);
				}
			}

			for ( auto key : toRemove ) {
				pmf.data.erase(key);
				cdf.data.erase(key);
			}
		}

		/**
			Quantises the distribution to fit a designated number of equi-spaced sample points.
		 */
		void Interpolate(int sample_points) {
			Histogram<double> h;

			double min = cdf.data.begin()->first ;
			double max = cdf.data.rbegin()->first;

			for ( int i = 0; i < sample_points; i++ ) {
				double t = min + i * (max - min) / ( sample_points - 1 );
				h.data[t] = 1000 * Cdf(t);
			}

			SetCdf( h );
		}
	};
}
