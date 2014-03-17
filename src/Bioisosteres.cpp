#include "Bioisosteres.hpp"
#include <iostream>
#include <vector>


namespace Bioisosteres{
	void testCall()
	{
		std::cout << "called function" << std::endl;
	}

	/*std::vector<sp_fragments> fragment_mol(ROMol& mol)
	{
		// take a dereferenced mol pointer and fragment
		ROMol* frag = MolFragmenter::fragmentOnBRICSBonds(mol);
		std::vector<sp_fragments> num_frags = MolOps::getMolFrags(*frag);
		return num_frags;
	}*/
}
