#pragma once
#include <functional>

#include "Exception.hpp"
#include "Function.hpp"

namespace Lmvq {
	class AbstractOptimiser {

	protected:
		Function & objectiveFunction;
		double * direction = 0;
		double * optimalX = 0;
		double * previousX = 0;
		size_t dim;
		double optimalY;
		double previousY;
		double epsilon = 0.001;
		uint maxIterations = 1000;
		uint iteration;
		std::function<void(void)> monitor;

		/// <summary> Create internal objects required to execute the
		/// search.
		/// </summary>

		void Initialise() {
			dim = objectiveFunction.InputDimension();
			optimalX = new double[dim];
			previousX = new double[dim];
			direction = new double[dim];
			Reset();
		}

	public:

		AbstractOptimiser(Function & objectiveFunction) : objectiveFunction(objectiveFunction) {
			Initialise();
		}

		virtual ~AbstractOptimiser() {
			delete[] optimalX;
			delete[] previousX;
			delete[] direction;
		}

		/// <summary> Get the objective function to be minimised.
		/// </summary>
		Function & ObjectiveFunction() {
			return objectiveFunction;
		}


		/// <summary> Set the initial parameter setting, updating xOpt, xPrev and direction
		/// vectors.
		/// </summary>

		void SetInitialX(double * value) {
			DAX::Set(optimalX, dim, value);
			Reset();
		}

		/// <summary> The most recently evaluated gradient of the objective function.
		/// </summary>

		double * Direction() { return direction; }

		/// <summary>
		/// The most recently evaluated estimate of the optimal input value.
		/// </summary>
		double * OptimalX() { return optimalX; }

		/// <summary>
		/// The most recent estimate of the optimal output value.
		/// </summary>
		double & OptimalY() { return optimalY; }

		/// <summary>
		/// The previously accepted input value.
		/// </summary>
		double * XPrev() { return previousX; }

		/// <summary>
		/// The previously accepted output value.
		/// </summary>
		double & YPrev() { return previousY; }

		/// <summary> Proximity-based convergence threshold:
		/// Converged if |xOpt - xPrev| < epsilon.
		/// </summary>

		double Epsilon() {
			return epsilon;
		}

		/// <summary> Proximity-based convergence threshold:
		/// Converged if |xOpt - xPrev| < epsilon.
		/// </summary>
		void SetEpsilon(double value) {
			if (value <= 0) {
				throw new Exception("", FileAndLine);
			}

			epsilon = value;
		}

		/// <summary>
		/// Returns the current iteration. This value is only meaningful 
		/// during execution of Run(), when it can be accessed by a Monitor 
		/// callback to provide feedback about optimisation progress.
		/// </summary>
		int Iteration() {
			return iteration;
		}

		/// <summary>
		/// The maximum number of steps to take in a given run.
		/// </summary>
		uint & MaxIterations() { return maxIterations; }

		/// <summary>
		/// Converged if |xOpt - xPrev| < epsilon.
		/// </summary>
		virtual bool Converged() {
			return DAX::DistanceSquared(optimalX, dim, previousX) < epsilon * epsilon;
		}

		/// <summary>
		/// A callback that lets us observe (and if necessary interfere with) the progress of
		///	a run.
		/// </summary>
		function<void(void)> & Monitor() {
			return monitor;
		}

		/// <summary>
		/// Restore the state to a condition from which it will try to 
		/// execute at least 1 step before testing for convergence.
		/// </summary>
		virtual void Reset() {
#pragma omp parallel for
			for ( size_t i = 0; i < dim; i++) {
				direction[i] = 0;
				previousX[i] = optimalX[i];
			}
			previousX[0] += 2 * epsilon;
			optimalY = numeric_limits<double>::max();
			previousY = numeric_limits<double>::max();
		}

		/// <summary> Steps maxIterations times, or until convergence, calling monitor (if) set
		/// after each iteration.
		/// <para>
		///		You will need to override Step and UpdateStepsize to 
		/// </para>
		/// </summary>

		void Run() {
			optimalY = objectiveFunction.Eval(optimalX);

			for (iteration = 1; iteration <= maxIterations; iteration++) {
				Step();

				optimalY = objectiveFunction.Eval(optimalX);

				if (Converged()) {
					break;
				}

				PostStepUpdate();

				if (monitor) {
					monitor();
				}
			}
		}

		/// <summary> Update xOpt according to the learning rule associated 
		/// with this optimiser.
		/// </summary>
		virtual void Step() = 0;

		/// <summary> Override this to update parameters such as step-sizes associated with 
		/// each component of the gradient vector. 
		/// </summary>

		virtual void PostStepUpdate() = 0;
	};
}
