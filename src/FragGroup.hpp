#ifndef FRAGGROUP_DEF

#define FRAGGROUP_DEF

// standard libraries
#include <iostream>
#include <vector>

// RDKit libraries
#include <GraphMol/ROMol.h>

// my libraries
#include "Frag.hpp"

class FragGroup
{
	public:
		// frag group
		std::vector<Frag> group;

		// constructor
		FragGroup(Frag frag);
		FragGroup(std::vector<Frag> new_group);
		
		// methods
		void addFragment(Frag fragment);
		void mergeFragGroup(FragGroup anotherGroup);

		int getLength();

	private:
		void getFragments();

};
#endif
