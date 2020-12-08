#pragma once

#include <FL/Fl.H>
#include <FL/Fl_Group.H>

namespace LBFL {
	class GroupMemory {
	private:
		const Fl_Group * current;

	public:
		GroupMemory() :
			current(Fl_Group::current()) //
		{}

		~GroupMemory() {
			Fl_Group::current(current);
		}
	};
}
