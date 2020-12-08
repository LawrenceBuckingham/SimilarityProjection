#pragma once

// Trick Visual studio.
#if __cplusplus < 201103L
#undef __cplusplus
#define __cplusplus 201703L
#endif

#include "EnumBase.hpp"

#include <string>
#include <mutex>
#include <vector>

using namespace std;

namespace QutBio {
	using Distance = int;

	typedef class DistanceType * pDistanceType;

	class DistanceType : public EnumBase {
	private:
		DistanceType(string literal, int value) : EnumBase(literal, value) {}
		static std::mutex m;

	public:
		static DistanceType * HalperinEtAl() {
			std::unique_lock < mutex > lck{ m };
			static DistanceType value("HalperinEtAl", 0);
			return &value;
		}

		static DistanceType * UngappedEdit() {
			std::unique_lock < mutex > lck{ m };
			static DistanceType value("UngappedEdit", 1);
			return &value;
		}

		static DistanceType * BlosumDistance() {
			std::unique_lock < mutex > lck{ m };
			static DistanceType value("BlosumDistance", 2);
			return &value;
		}

		/**
		 *	<summary>
		 *		A custom distance which is defined by a similarity matrix loaded at run time.
		 *	</summary>
		 */
		static DistanceType * Custom() {
			std::unique_lock < mutex > lck{ m };
			static DistanceType value("Custom", 3);
			return &value;
		}

		static vector<EnumBase *> Values() {
			std::vector<EnumBase *> result(4);
			result[0] = HalperinEtAl();
			result[1] = UngappedEdit();
			result[2] = BlosumDistance();
			result[3] = Custom();
			return result;
		}
	};

	typedef DistanceType *pDistanceType;
}
