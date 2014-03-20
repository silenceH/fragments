#ifndef BIOISOSTERES_DEF
#define BIOISOSTERES_DEF
/* 
 * class for generating bioisosteres
 */

// standard libraries
#include <iostream>

// RDKit libraries
#include <GraphMol/ROMol.h>

using namespace RDKit;

typedef boost::shared_ptr<ROMol> sp_fragments;

namespace Bioisosteres
{
	void testCall();
	std::vector<sp_fragments> fragment_mol(ROMol& mol);
	std::vector<ROMol*> getMolsFromFile(std::string file_name);

}
#endif
