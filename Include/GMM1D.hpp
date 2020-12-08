#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <vector>
#include "Array.hpp"
#include "Assert.hpp"
#include "Distribution.hpp"
#include "NormalDistribution.hpp"
#include "Types.hpp"
#include <omp.h>

namespace QutBio {
	/***
	 *	Statistical information from:
	 *	http://stats.stackexchange.com/questions/16608/what-is-the-variance-of-the-weighted-mixture-of-two-gaussians
	 *
	 *	https://github.com/szpiech/gmm
	 *
	 *	http://ssli.ee.washington.edu/people/bilmes/mypapers/em.pdf
	 ***/
	class GMM1D : public  Distribution {
	private:
		std::vector<double> a;
		std::vector<double> mu;
		std::vector<double> sigma;
		double aicc = 0;
		double epsilon;

	public:
		/**
		**	<summary>
		**		Initialises the mixture of Gaussians to hold the specified number of
		**		kernels. Each is placed at a random, small displacement from zero, with
		**		bandwidth of 1.0.
		**	</summary>
		**	<param name="epsilon">Bound used to terminate the bisection search.
		**	</param>
		*/
		GMM1D(size_t size = 1, double epsilon = 1e-5) : epsilon(epsilon) {
			for ( size_t i = 0; i < size; i++ ) {
				mu.push_back(1e-5 * (double) rand() / RAND_MAX);
				sigma.push_back(1);
				a.push_back(1.0 / size);
			}
		}

		/***
		**	<summary>
		**		Initialises the mixture of Gaussians with kernel parameters provided in a
		**		collection of parallel arrays.
		**	</summary>
		**	<param name="a">
		**		Vector of mixing weights. These should be normalised so that
		**		they sum to 1.
		**	</param>
		**	<param name="mu">
		**		Vector containing centres of kernels.
		**	</param>
		**	<param name="sigma">
		**		Vector containing bandwidths of kernels.
		**	</param>
		**	<param name="epsilon">
		**		Bound used to terminate the bisection search in InverseCdf.
		**	</param>
		*/
		GMM1D(
			const std::vector<double> & a,
			const std::vector<double> & mu,
			const std::vector<double> & sigma,
			double aicc,
			double epsilon = 1e-5
			) : a(a), mu(mu), sigma(sigma), aicc(aicc), epsilon(epsilon) {
		}

		/***
		 * Gets the squared Euclidean distance from the parameters of this
		 * instance to those of another, which is required to have the same
		 * internal dimensions as this.
		 ***/
		double Distance(const GMM1D & other) {
			double d2 = 0;

			assert_intsEqual(mu.size(), other.mu.size());

			for ( uint i = 0; i < mu.size(); i++ ) {
				double t = a[i] - other.a[i]; d2 += t * t;
				t = mu[i] - other.mu[i]; d2 += t * t;
				t = sigma[i] - other.sigma[i]; d2 += t * t;
			}

			return d2;
		}

		/***
		 * Gets the squared Euclidean distance from the parameters of this
		 * instance to those of another, which is required to have the same
		 * internal dimensions as this.
		 ***/
		double Distance(
			const std::vector<double> & other_a,
			const std::vector<double> &other_mu,
			const std::vector<double> & other_sigma
			) {
			double d2 = 0;

			assert_intsEqual(mu.size(), other_mu.size());

			for ( uint i = 0; i < mu.size(); i++ ) {
				double t = a[i] - other_a[i];
				d2 += t * t;
				t = mu[i] - other_mu[i];
				d2 += t * t;
				t = sigma[i] - other_sigma[i];
				d2 += t * t;
			}

			return d2;
		}

		template<typename T>
		void Initialise(const std::vector<T> & data) {
			for ( uint i = 0; i < a.size(); i++ ) {
				double t = data[rand() % data.size()];
				mu[i] = t;
			}
		}

		/***
		 *	<summary>
		 *		Apply the EM algorithm to fit a mixture model to a set of numeric data.
		 *	</summary>
		 *	<param name="data">
		 *		A vector of numeric values, that is, they can be converted automatically
		 *		to double by compiler generated code.
		 *	</param>
		 *	<param name="epochs">
		 *		The number of times to present each training observation.
		 *	</param>
		 *	<param name="epsilon">
		 *		A positive floating point termination threshold. training will cease when
		 *		the squared Euclidean distance between successive iterations of the model's parameter
		 *		vectors shrinks below the value.
		 *	</param>
		 *	<typeparam name="T">
		 *		The numeric type of the data to be modelled by this mixture.
		 *	</typeparam>
		 ***/

		template<typename T>
		void Train(const std::vector<T> & data, const uint epochs, const double epsilon, bool verbose = false) {
			using std::vector;
			size_t m = a.size();
			size_t n = data.size();
			const T *x = data.data();

			vector<double> a2(a);
			vector<double> mu2(mu);
			vector<double> sigma2(sigma);

			double * curr_a(a.data());
			double * curr_mu(mu.data());
			double * curr_sigma(sigma.data());

			double * next_a(a2.data());
			double * next_mu(mu2.data());
			double * next_sigma(sigma2.data());

			double ** p;
			double ** sum_p;
			double ** sum_p_x;
			double ** sum_p_x_squared;

			RawMatrix<double> pijm(n, m);
			vector<vector<double>> & pij( pijm.data );

			int numThreads = 0;
			time_t start_time = time(0);

#if USE_OMP
#pragma omp parallel
#endif
			numThreads = omp_get_num_threads();

			if ( verbose ) {
				cerr << "Training Gaussian mixture model with " << numThreads << " threads." << endl;
				cerr << "Start time: " << start_time << endl;
			}

			p = (double **) alloca(m*sizeof(double *));
			sum_p = (double **) alloca(m*sizeof(double *));
			sum_p_x = (double **) alloca(m*sizeof(double *));
			sum_p_x_squared = (double **) alloca(m*sizeof(double *));

			// Create accumulators for each thread.
			for ( size_t j = 0; j < m; j++ ) {
				p[j] = (double *) alloca(numThreads*sizeof(double));
				sum_p[j] = (double *) alloca(numThreads*sizeof(double));
				sum_p_x[j] = (double *) alloca(numThreads*sizeof(double));
				sum_p_x_squared[j] = (double *) alloca(numThreads*sizeof(double));
			}

			auto swap_pointers = [&](void) {
				double *t;
				t = next_a; next_a = curr_a; curr_a = t;
				t = next_mu; next_mu = curr_mu; curr_mu = t;
				t = next_sigma; next_sigma = curr_sigma; curr_sigma = t;
			};

			for ( uint epoch = 0; epoch < epochs; epoch++ ) {
				for ( uint j = 0; j < m; j++ ) {
					for ( int thread = 0; thread < numThreads; thread++ ) {
						// Clear the accumulators for each parameter.
						p[j][thread] = 0;
						sum_p[j][thread] = 0;
						sum_p_x[j][thread] = 0;
						sum_p_x_squared[j][thread] = 0;
					}
				}

				// calculate sum of probabilities, a and mu for next step.
#if USE_OMP
#pragma omp parallel for
#endif
				for (size_t i = 0; i < n; i++ ) {
					int thread = omp_get_thread_num();

					// This is local to each training observation (indexed by distinct i values)
					// so it is not subject to race conditions.
					double total_p = 0;

					for ( uint j = 0; j < m; j++ ) {
						double t = max(1e-10, NormalDistribution::Pdf(x[i], curr_mu[j], curr_sigma[j]) * curr_a[j]);
						pij[i][j] = t;

						// TODO - Do something more elegant to permit the probabilities to be evaluated
						// TODO - Similar to what I did in LMVQ.
						if ( !isfinite(t) ) {
							p[j][thread] = 0;
							continue;
						}

						p[j][thread] = t;
						total_p += t;
					}

					for ( uint j = 0; j < m; j++ ) {
						double w = p[j][thread] / total_p;
						sum_p[j][thread] += w;
						sum_p_x[j][thread] += x[i] * w;
					}
				}

				for ( uint j = 0; j < m; j++ ) {
					// Collect the partial sums computed by individual threads into sum.
					double sum_p_j = 0, sum_p_x_j = 0;

					for ( int thread = 0; thread < numThreads; thread++ ) {
						sum_p_j += sum_p[j][thread];
						sum_p_x_j += sum_p_x[j][thread];
					}

					next_a[j] = sum_p_j / n;
					next_mu[j] = sum_p_x_j / sum_p_j;
				}

				// calculate variance.
				// TODO - Optimise this if it needs optimising.
#if USE_OMP
#pragma omp parallel for
#endif
				for (size_t i = 0; i < n; i++ ) {
					int thread = omp_get_thread_num();

					double total_p = 0;

					for ( uint j = 0; j < m; j++ ) {
						double t = pij[i][j];

						// TODO - Do something more elegant to permit the probabilities to be evaluated
						// TODO - Similar to what I did in LMVQ.
						if ( !isfinite(t) ) continue;

						p[j][thread] = t;
						total_p += t;
					}

					if ( total_p > 0 ) {

						for ( uint j = 0; j < m; j++ ) {
							double t = x[i] - next_mu[j];
							sum_p_x_squared[j][thread] += t * t * p[j][thread] / total_p;
						}
					}
				}

				for ( uint j = 0; j < m; j++ ) {
					double sum_p_j = 0;
					double sum_p_x_squared_j = 0;

					for ( int thread = 0; thread < numThreads; thread++ ) {
						sum_p_j += sum_p[j][thread];
						sum_p_x_squared_j += sum_p_x_squared[j][thread];
					}

					next_sigma[j] = sqrt(sum_p_x_squared_j / sum_p_j);
					//cerr << "j = " << j << ", a[j] = " << next_a[j] 
					//	<< ", mu[j] = " << next_mu[j]
					//	<< ", sigma[j] = " << next_sigma[j]
					//	<< endl;
				}

				if ( verbose ) {
					cerr << "Epoch " << epoch << ": LogLikelihood = " << LogLikelihood(data) << endl;
				}

				if ( Distance(a2, mu2, sigma2) < epsilon ) {
					break;
				}

				swap_pointers();
			}

			if ( curr_a == a.data() ) {
				a = a2;
				mu = mu2;
				sigma = sigma2;
			}

			if ( verbose ) {
				time_t end_time = time(0);
				cerr << "End time: " << end_time << endl;
				cerr << "Elapsed time: " << (end_time - start_time) << endl;
			}
		}

		double Cdf(double t) {
			double total_p = 0;
			const size_t m = a.size();

			for ( uint j = 0; j < m; j++ ) {
				double pj = NormalDistribution::Cdf(t, mu[j], sigma[j]) * a[j];

				// TODO - Do something more elegant to permit the probabilities to be evaluated
				// TODO - Similar to what I did in LMVQ.
				if ( !isfinite(pj) ) {
					continue;
				}

				total_p += pj;
			}

			return total_p;
		}

		double Pdf(double t) {
			double total_p = 0;
			const size_t m = a.size();

			for ( uint j = 0; j < m; j++ ) {
				double pj = NormalDistribution::Pdf(t, mu[j], sigma[j]) * a[j];

				// TODO - Do something more elegant to permit the probabilities to be evaluated
				// TODO - Similar to what I did in LMVQ.
				if ( !isfinite(pj) ) {
					continue;
				}

				total_p += pj;
			}

			return total_p;
		}

		/***
		*	<summary>
		*		Gets the natural logarithm of the likelihood of a set of numeric values
		*		given the current parameters of the model.
		*	</summary>
		*	<param name="sample">
		*		A vector of numeric values, that is, they can be converted automatically
		*		to double by compiler generated code.
		*	</param>
		*	<typeparam name="T">
		*		The numeric type of the data to be modelled by this mixture.
		*	</typeparam>
		***/

		template<typename T>
		double LogLikelihood(const std::vector<T> & sample) {
			int numThreads;

#pragma omp parallel
			numThreads = omp_get_num_threads();

			double logLikelihood = 0;
			double * ll_per_thread = (double *) alloca(numThreads * sizeof(double));

			for ( int t = 0; t < numThreads; t++ ) {
				ll_per_thread[t] = 0;
			}

#pragma omp parallel for
			for (size_t i = 0; i < sample.size(); i++ ) {
				int t = omp_get_thread_num();
				ll_per_thread[t] += log(Pdf(sample[i]));
			}

			for ( int t = 0; t < numThreads; t++ ) {
				logLikelihood += ll_per_thread[t];
			}

			return logLikelihood;
		}

		/***
		*	<summary>
		*		Gets the corrected Akaike Information Criterion of the model for a set of numeric
		*		values.
		*		The value is also stored in the model for later reference.
		*	</summary>
		*	<param name="sample">
		*		A vector of numeric values, that is, they can be converted automatically
		*		to double by compiler generated code.
		*	</param>
		*	<typeparam name="T">
		*		The numeric type of the data to be modelled by this mixture.
		*	</typeparam>
		***/

		template<typename T>
		double AICc(const std::vector<T> & sample) {
			size_t k = 3 * mu.size(); // a, mu, sigma -> three parameters for each component.
			size_t n = sample.size();
			double AIC = 2 * k - 2 * LogLikelihood(sample);
			return aicc = (n > (k + 1) ? AIC + 2 * k * (k + 1) / (n - k - 1) : numeric_limits<double>::max());
		}

		/**
		 *	Gets the Akaike Information Criterion stored in the model.
		 */
		double AICc() { return aicc; }

		/**
		 *	Gets the Akaike Information Criterion stored in the model.
		 */
		void SetAICc( double aicc ) { this->aicc = aicc; }

		/**
		**	<summary>
		**		Brute-force implementation of inverse CDF by first bracketing
		**		the inverse value in an interval, and subsequently using 
		**		bisection search to locate it.
		**	</summary>
		**	<param name="t">The probability threshold at which the inverse value 
		**		is desired.
		**	</param>
		**	<returns>
		**		Value p such that Cdf(p-epsilon/2) <= t and Cdf(p+epsilon/2) >= t.
		**	</returns>
		*/
		double InverseCdf(double t) {
			if ( t <= 0 ) return numeric_limits<double>::lowest();
			else if ( t >= 1 ) return numeric_limits<double>::max();
			else {
				// Start with mean as initial bracket endpoint
				double a = Mean();
				double fa = Cdf(a);

				// choose another endpoint 1 unit in the appropriate direction.
				double dir = t < fa ? -1 : +1;

				double b = a + dir;
				double fb = Cdf(b);

				double lo = min(a, b);
				double hi = max(a, b);

				double fLo = min(fa, fb);
				double fHi = max(fa, fb);

				while ( t < fLo ) {
					lo += lo - hi;
					fLo = Cdf(lo);
				}

				while ( t > fHi ) {
					hi += hi - lo;
					fHi = Cdf(hi);
				}

				double mid = 0;

				while ( hi - lo > 1e-5 ) {
					mid = (lo + hi) / 2;
					double fMid = Cdf(mid);

					if ( t < fMid ) {
						hi = mid;
						fHi = fMid;
					}
					else {
						lo = mid;
						fLo = fMid;
					}
				}

				return mid;
			}
		}

		double Mean(void) {
			double m = 0;

			for ( uint i = 0; i < a.size(); i++ ) {
				m += a[i] * mu[i];
			}

			return m;
		}

		double StdDev(void) {
			double variance = 0;

			for ( uint i = 0; i < a.size(); i++ ) {
				variance += a[i] * (sigma[i] * sigma[i] + mu[i] * mu[i]);
			}

			double m = Mean();

			variance -= m * m;

			return sqrt(variance);
		}

		size_t Size() { return mu.size(); }

		friend std::ostream & operator<<(std::ostream & out, GMM1D & model) {
			out << "alpha,mu,sigma" << endl;

			for ( size_t i = 0; i < model.mu.size(); i++ ) {
				out << model.a[i] << "," << model.mu[i] << "," << model.sigma[i] << endl;
			}

			out << endl;
			return out;
		}

		static void Parse(std::istream & in, std::vector<GMM1D> & result ) {
			using std::string;
			using std::getline;
			using std::vector;

			string s;

			while ( std::getline(in, s) ) {
				double aicc = 0;

				if ( s.find( "AIC" ) == 0 )
				{
					string ignored;
					istringstream str(s);
					str >> ignored >> aicc;
				}
				if ( s.find("alpha,mu,sigma") != s.npos ) {
					vector<double> alpha;
					vector<double> mu;
					vector<double> sigma;

					while ( std::getline(in, s) && s.length() > 0 ) {
						double alpha_, mu_, sigma_;
						char c1, c2;
						istringstream str(s);
						str >> alpha_ >> c1 >> mu_ >> c2 >> sigma_;
						alpha.push_back(alpha_);
						mu.push_back(mu_);
						sigma.push_back(sigma_);
					}

					GMM1D model( alpha, mu, sigma, aicc );
					result.push_back( model );
				}
			}
		}
	};
}
