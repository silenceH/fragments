// standard libraries
#include <iostream> 

// RDKit libraries 
#include <GraphMol/ROMol.h>

// FragGroup header
#include "FragGroup.hpp"
#include "Frag.hpp"


using namespace RDKit; 

// constructors
FragGroup::FragGroup(Frag frag)
{
	group.push_back(frag);
}

FragGroup::FragGroup(std::vector<Frag> new_group)
{
	group = new_group;
}

// methods
int FragGroup::getLength()
{
	return group.size();
}
