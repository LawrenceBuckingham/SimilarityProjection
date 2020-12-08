#pragma once 
#include "Types.hpp"
#include "Mapping.hpp"

// TODO: Copyright, license and documentation.
namespace Lmvq {
	class Function : public Mapping {
	public:
		virtual void Map(double x[], double y[]) {
			y[0] = Eval(x);
		}

		virtual uint OutputDimension(void) {
			return 1;
		}

		virtual double Eval(double x[]) = 0;
		virtual void GetGradient(double x[], double grad[]) = 0;
	};
}
