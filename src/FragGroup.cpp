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

FragGroup::FragGroup(std::vector<sp_fragments> new_group)
{
	for (const auto& fr : new_group)
	{
		group.push_back(Frag(fr));
	}
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
