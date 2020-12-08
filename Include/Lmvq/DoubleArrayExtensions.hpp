#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include <iostream> 
#include <string> 
#include <sstream> 
#include <functional>
#include <cstring>
#include <cmath>

namespace Lmvq {

	/// <summary>
	///		The extension methods provided here permit raw arrays of double to 
	///		be treated like AbstractVector objects, but with reduced overhead.
	/// </summary>

	class DoubleArrayExtensions {

	public:

		/// <summary>
		///		Replaces the contents of the destination vector with that of 
		///		the source vector.
		/// </summary>
		/// <param name="destination">
		///		The array into which values are to be copied.
		/// </param>
		/// <param name="source">
		///		the array from which values will be copied.
		/// </param>
		/// <returns>
		///		Returns a copy of the destination array, to permit convenient 
		///		chaining of method calls.
		/// </returns>

		static double * Set(
			double * destination,
			size_t n,
			const double * source
			) {
			std::memcpy(destination, source, n * sizeof(double));
			return destination;
		}

		/// <summary>
		///		Replaces the contents of the destination vector with the scalar 
		///		value provided.
		/// </summary>
		/// <param name="destination">
		///		The array into which values are to be copied.
		/// </param>
		/// <param name="scalar">
		///		The scalar value that is to fill the array.
		/// </param>
		/// <returns>
		///		Returns a copy of the destination array, to permit convenient 
		///		chaining of method calls.
		/// </returns>

		static double * Set(
			double * destination,
			size_t n,
			double scalar
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] = scalar;
			}

			return destination;
		}

		/// <summary> Duplicates the supplied source array.
		/// </summary>
		/// <param name="source">
		///		An array that is to be duplicated.
		/// </param>
		/// <returns>
		///		A duplicate copy of the source array.
		/// </returns>

		static double * Duplicate(double * source, int n) {
			double * result = new double[n];
			Set(result, n, source);
			return result;
		}

		/// <summary> Replaces the contents of this vector with a scaled set of 
		///		values from a list.
		/// </summary>
		/// <param name="a">The scalar multiplier.</param>
		/// <param name="result">The list of values.</param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Set(
			double * destination,
			size_t n,
			double a,
			double * x
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] = a * x[i];
			}

			return destination;
		}

		/// <summary> Replaces the contents of this vector with a weighted set 
		///		of values drawn from a list.
		/// </summary>
		/// <param name="a">The list of weights.</param>
		/// <param name="result">The list of values.</param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Set(
			double * destination,
			size_t n,
			double * a,
			double * x
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] = a[i] * x[i];
			}

			return destination;
		}

		/// <summary> Adds a vector to this vector.
		/// <para>That is, this += result</para>
		/// </summary>
		/// <param name="other">The direction vector.</param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Add(
			double * destination,
			size_t n,
			double * other
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] += other[i];
			}

			return destination;
		}

		/// <summary> Adds a scaled vector from this vector.
		/// <para>That is, this += a * result</para>
		/// </summary>
		/// <param name="a">The scalar multiplier.</param>
		/// <param name="result">The direction vector.</param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Add(
			double * destination,
			size_t n,
			double a,
			double * x
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] += a * x[i];
			}

			return destination;
		}

		/// <summary> Adds a scaled vector from this vector.
		/// <para>That is, this += a * result</para>
		/// </summary>
		/// <param name="a">The scalar multiplier.</param>
		/// <param name="result">The direction vector.</param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		template<class T>
		static double * Add(
			double * destination,
			size_t n,
			double a,
			const T & x
			) {
			int i = 0;

			for ( double t : x ) {
				if ( i >= n ) break;
				destination[i++] += a * t;
			}

			return destination;
		}

		/// <summary> Subtracts a vector from this vector.
		/// <para>That is, this -= result</para>
		/// </summary>
		/// <param name="other">The direction vector.</param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Sub(
			double * destination,
			size_t n,
			double * other
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] -= other[i];
			}

			return destination;
		}

		/// <summary> Subtracts a scaled vector from this vector.
		/// <para>That is, this -= a * result</para>
		/// </summary>
		/// <param name="a">The scalar multiplier.</param>
		/// <param name="result">The direction vector.</param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Sub(
			double * destination,
			size_t n,
			double a,
			double * x
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] -= a * x[i];
			}

			return destination;
		}

		/// <summary> Negates the destination vector, element by element.
		/// </summary>
		/// <param name="destination">
		///		The array to be negated.
		/// </param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Negate(
			double * destination,
			size_t n
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] = -destination[i];
			}

			return destination;
		}

		/// <summary> Computes the scalar product of the subject vector with another. 
		/// </summary>
		/// <param name="subject">
		///		The subject vector.
		/// </param>
		/// <param name="other">
		///		The other vector.
		/// </param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double Dot(
			double * subject,
			size_t n,
			double * other
			) {
			double result = 0;

			for ( size_t i = 0; i < n; i++ ) {
				result += subject[i] * other[i];
			}

			return result;
		}

		/// <summary> Multiplies each element of the destination array by a scalar.
		/// </summary>
		/// <param name="destination">
		///		The destination array.
		/// </param>
		/// <param name="scalar">
		///		A scalar multiplier.
		/// </param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Mul(
			double * destination,
			size_t n,
			double scalar
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] *= scalar;
			}

			return destination;
		}

		/// <summary> Multiplies each element of the destination array by the
		///		corresponding element of an array of weights.
		/// </summary>
		/// <param name="destination">
		///		The destination array.
		/// </param>
		/// <param name="weights">
		///		An array of element weights which will be used to scale the
		///		corresponding elements of the destination array.
		/// </param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Mul(
			double * destination,
			size_t n,
			double * weights
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] *= weights[i];
			}

			return destination;
		}

		/// <summary> Divides each element of the destination array by a scalar.
		/// </summary>
		/// <param name="destination">
		///		The destination array.
		/// </param>
		/// <param name="scalar">
		///		A scalar divisor.
		/// </param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Div(
			double * destination,
			size_t n,
			double scalar
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] /= scalar;
			}

			return destination;
		}

		/// <summary> Divides each element of the destination array by the
		///		corresponding element of an array of weights.
		/// </summary>
		/// <param name="destination">
		///		The destination array.
		/// </param>
		/// <param name="weights">
		///		An array of element weights which will be used to scale the
		///		corresponding elements of the destination array.
		/// </param>
		/// <returns>
		///		The destination array for convenient composition of operations.
		///	</returns>

		static double * Div(
			double * destination,
			size_t n,
			double * weights
			) {
			for ( size_t i = 0; i < n; i++ ) {
				destination[i] /= weights[i];
			}

			return destination;
		}

		/// <summary> Test to see if the subject array is more similar (in terms of 
		///		squared Euclidean distance) to another than some previous item.
		/// </summary>
		/// <param name="other">
		///		The other array.
		/// </param>
		/// <param name="subject">
		///		The subject array.
		/// </param>
		/// <param name="previousDistSquared">
		///		The previous best (smallest) squared Euclidean distance. 
		/// </param>
		/// <returns>
		///		True iff ||subject-other||^2 &lt; previousDistanceSquared_0, and
		///		also if the result is true, previousdistanceSquared_1 = 
		///		||subject-other||^2.
		/// </returns>

		static bool IsNearer(
			double * subject,
			size_t n,
			double * other,
			double & previousDistSquared
			) {
			double distanceSquared = 0;

			for ( size_t i = 0; i < n; i++ ) {
				double t = subject[i] - other[i];

				distanceSquared += t * t;

				if ( distanceSquared >= previousDistSquared ) {
					return false;
				}
			}

			previousDistSquared = distanceSquared;

			return true;
		}

		static double Distance(
			double * subject,
			size_t n,
			double * other
			) {
			return sqrt(DistanceSquared(subject, n, other));
		}

		static double DistanceSquared(
			double * subject,
			size_t n,
			double * other
			) {
			double distanceSquared = 0;

			for ( size_t i = 0; i < n; i++ ) {
				double t = subject[i] - other[i];

				distanceSquared += t * t;
			}

			return distanceSquared;
		}

		/// <summary> Compute the (squared?) Mahalanobis distance to a point.
		/// </summary>
		/// <param name="other"></param>
		/// <param name="inverseCovariance"></param>
		/// <returns></returns>

		static double Distance(
			double * subject,
			size_t n,
			double * other,
			double ** inverseCovariance
			) {
			double sum = 0;

			for ( size_t i = 0; i < n; i++ ) {
				double t = subject[i] - other[i];

				for ( size_t j = 0; j < n; j++ ) {
					sum += inverseCovariance[i][j] * t * (subject[j] - other[j]);
				}
			}

			return sum;
		}

		// Not totally satisfactory, but it'll do for now.
		static bool IsEqualTo(
			double * subject,
			size_t n,
			double * other
			) {
			if ( subject == other ) {
				return true;
			}

			for ( size_t i = 0; i < n; i++ ) {
				if ( subject[i] != other[i] ) {
					return false;
				}
			}

			return true;
		}

		static double NormSquared(
			double * subject,
			size_t n
			) {
			return Dot(subject, n, subject);
		}

		static double Norm(
			double * subject,
			size_t n
			) {
			return std::sqrt(NormSquared(subject, n));
		}

		static double * Apply(
			double * subject,
			size_t n,
			std::function<void(double)> f
			) {
			for (size_t i = 0; i < n; i++ ) {
				f(subject[i]);
			}

			return subject;
		}

		/// <summary> Applies the supplied function to each element of the subject 
		///		array, replacing the elements with the resulting values.
		/// </summary>
		/// <param name="subject"></param>
		/// <param name="f"></param>
		/// <returns>
		///		The subject array.
		/// </returns>

		static double * Update(
			double * subject,
			size_t n,
			std::function<double(double)> f
			) {
			for (size_t i = 0; i < n; i++ ) {
				subject[i] = f(subject[i]);
			}

			return subject;
		}

		static void HardLimit(
			double * subject,
			size_t n
			) {
			Update(subject, n, [&](double x) { return round(x); });
		}

		static void HardLimit(
			double ** subject,
			size_t m,
			size_t n
			) {
			for (size_t i = 0; i < m; i++ ) {
				HardLimit(subject[i], n);
			}
		}

		static void CopyTo(
			double * subject,
			size_t n,
			double * x
			) {
			for (size_t i = 0; i < n; i++ ) {
				x[i] = subject[i];
			}
		}

		/// <summary> Compares two vectors.
		/// </summary>
		/// <param name="other"></param>
		/// <returns></returns>

		static int CompareTo(
			double * subject,
			size_t n,
			double * other
			) {
			int result = 0;

			for (size_t i = 0; result == 0 && i < n; i++ ) {
				result = subject[i] < other[i] ? -1 : subject[i] > other[i] ? 1 : 0;
			}

			return result;
		}

		static string ConvertToString(
			double * subject,
			size_t n
			) {
			std::ostringstream w;
			w << '[';

			if ( n > 0 ) {
				w << subject[0];

				for (size_t i = 1; i < n; i++ ) {
					w << ',' << subject[i];
				}
			}

			w << ']';
			return w.str();
		}

	};
}
