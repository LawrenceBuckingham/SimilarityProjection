#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "AbstractOptimiser.hpp"
#include "../Exception.hpp"

namespace Lmvq {

	using DAX = DoubleArrayExtensions;

	class SilvaAlmeidaOptimiser : public AbstractOptimiser {
	protected:
		bool rollback = false;
		double * stepSize = 0;
		double * previousDirection = 0;
		double stepSizeMultiplier = 1.0;
		double stepSizeDivisor = 2.0;
		double initialStepSize = 0.1;

		/// <summary> Executes the base class initialise procedure, then 
		/// initialises stepSize and previousDirection arrays.
		/// <para>
		/// Precondition: the objective function must have been set.
		/// Post-condition: stepSize and previousDirection have been allocated.
		/// </para>
		/// </summary>

		void Initialise() {
			previousDirection = new double[dim];
			stepSize = new double[dim];
			DAX::Set(stepSize, dim, initialStepSize);
		}

	public:
		SilvaAlmeidaOptimiser(Function & objectiveFunction) : AbstractOptimiser(objectiveFunction) {
			Initialise();
		}

		virtual ~SilvaAlmeidaOptimiser(void) {
			if ( stepSize ) delete[] stepSize;
			if ( previousDirection ) delete[] previousDirection;
		}

		/// <summary> Get the initial step size. Best called before 
		/// the InitialX mutator is called.
		/// </summary>

		double InitialStepSize() { return initialStepSize; }

		/// <summary> Get the initial step size. Best called before 
		/// the InitialX mutator is called.
		/// </summary>
		/// <exception cref="ArgumentExcception">
		/// Thrown if value is less than or equal to zero.
		/// </exception>

		void SetInitialStepSize(double value) {
			if ( value <= 0 ) throw Exception("Initial step size must be greater than zero.", FileAndLine);
			initialStepSize = value;
		}


		/// <summary> Gets the previous direction vector.
		/// </summary>


		double * PreviousDirection() {
			return previousDirection;
		}

		/// <summary> Gets the roll-back status of the optimiser.
		/// </summary>

		bool Rollback() {
			return rollback;
		}

		/// <summary> The step size to be used. This is adjusted dynamically as
		/// the program runs.
		/// </summary>

		double * StepSize() {
			return stepSize;
		}

		/// <summary> Get or set the step size divider, used to
		/// decrease stepSizes whenever the objective function 
		/// increases or the sign of the partial derivative changes.
		/// </summary>
		/// <exception cref="ArgumentException">
		///		If the supplied value is less than or equal to 1.0.
		/// </exception>

		double StepSizeDivisor() {
			return stepSizeDivisor;
		}

		/// <summary> Get or set the step size divider, used to
		/// decrease stepSizes whenever the objective function 
		/// increases or the sign of the partial derivative changes.
		/// </summary>
		/// <exception cref="ArgumentException">
		///		If the supplied value is less than or equal to 1.0.
		/// </exception>

		void SetStepSizeDivisor(double value) {
			if ( value <= 1.0 ) {
				throw Exception("StepSizeDivisor must be >= 1.", FileAndLine);
			}

			stepSizeDivisor = value;
		}

		/// <summary> Get or set the step size multiplier, used to
		/// increase stepSizes whenever the objective function 
		/// decreases
		/// </summary>
		/// <exception cref="ArgumentException">
		///		If the supplied value is less than or equal to 1.0.
		/// </exception>

		double StepSizeMultiplier() {
			return stepSizeMultiplier;
		}

		/// <summary> Get or set the step size multiplier, used to
		/// increase stepSizes whenever the objective function 
		/// decreases
		/// </summary>
		/// <exception cref="ArgumentException">
		///		If the supplied value is less than or equal to 1.0.
		/// </exception>

		void SetStepSizeMultiplier(double value) {
			if ( value <= 1.0 ) {
				throw Exception("StepSizeMultiplier must be >= 1.", FileAndLine);
			}

			stepSizeMultiplier = value;
		}

		/// <summary> Restore the state to a condition from which it will try to 
		/// execute at least 1 step before testing for convergence.
		/// </summary>

		void Reset() {
			AbstractOptimiser::Reset();
			rollback = false;
		}

		/// <summary> Evaluates the gradient and takes a step in that direction 
		///	(positive if maximising, negative if minimising).
		/// </summary>

		void Step() {
			if ( !rollback ) {
#pragma omp parallel for
				for (size_t i = 0; i < dim; i++ ) {
					previousX[i] = optimalX[i];
					previousDirection[i] = direction[i];
				}
				previousY = optimalY;
				objectiveFunction.GetGradient(optimalX, direction);
			}

#pragma omp parallel for
			for (size_t i = 0; i < dim; i++ ) {
				optimalX[i] -= stepSize[i] * direction[i];
			}
		}

		/// <summary> Check for roll-back and update step sizes.
		/// </summary>

		void PostStepUpdate() {
			rollback = optimalY >= previousY;

			if ( rollback ) {
#pragma omp parallel for
				for (size_t i = 0; i < dim; i++ ) {
					stepSize[i] /= stepSizeMultiplier;
					optimalX[i] = previousX[i];
				}
				optimalY = previousY;
			}
			else {
#pragma omp parallel for
				for (size_t i = 0; i < dim; i++ ) {
					if ( previousDirection[i] * direction[i] >= 0 ) {
						stepSize[i] *= stepSizeMultiplier;
					}
					else {
						stepSize[i] /= stepSizeMultiplier;
					}
				}
			}
		}
	};
}
