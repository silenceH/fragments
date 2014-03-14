#include <vector>
#include <stdlib.h> // exit getenv
#include "RDKMols.hpp"

std::vector<sp_fragments> RDKMols::fragment_mol(ROMol& mol)
{
	// take a dereferenced mol pointer and fragment
	ROMol* frag = MolFragmenter::fragmentOnBRICSBonds(mol);
	std::vector<sp_fragments> num_frags = MolOps::getMolFrags(*frag);
	return num_frags;
}


/* 
 * Method that takes a file name and returns a std::vector of ROMol pointers
 */
std::vector<ROMol*> RDKMols::getMols(std::string file_name)
{
	char* data = getenv("DATA");
	if (data==NULL)
	{
		std::cerr << "cannot find data environment" << std::endl;
		exit(EXIT_FAILURE);
	}
	

	std::string fname = std::string(data) + std::string("validation_overlays/") +file_name + std::string(".sdf");
	SDMolSupplier suppl(fname,true);	// sanitize mols and keep H atoms
	std::vector<ROMol*> mols;
	while(!suppl.atEnd())
	{
		ROMol *nmol = suppl.next();
		if(nmol!=0) 			// NULL is 0 in c++
			mols.push_back(nmol);
	}
	return mols;
}
