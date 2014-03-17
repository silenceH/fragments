#include <iostream>
#include "Frag.hpp"
#include "Bioisosteres.hpp"
#include "Bioisosteres.cpp"
#include "RDKMols.hpp"
#include "RDKMols.cpp"


int main() 
{
	Frag* new_frag1 = new Frag;
	Frag* new_frag2 = new Frag;
	
	new_frag1->setSmiles("smiles1");
	new_frag2->setSmiles("smiles2");
	std::string smiles1 = new_frag1->getSmiles();
	std::string smiles2 = new_frag2->getSmiles();
	std::cout << smiles1 << std::endl;
	std::cout << smiles2 << std::endl;

	Frag new_frag3;
	Frag new_frag4;
	
	new_frag3.setSmiles("smiles1");
	new_frag4.setSmiles("smiles2");
	std::string smiles3 = new_frag3.getSmiles();
	std::string smiles4 = new_frag4.getSmiles();
	std::cout << smiles3 << std::endl;
	std::cout << smiles4 << std::endl;

	// get mols 
	RDKMols* file = new RDKMols(std::string("P39900"));

	std::cout << (*file).ligand << std::endl;
	
	std::cout << "Number of molecules: " << file->mols.size() << std::endl;
	

	for(const auto& mol : file->mols)
	{
		std::cout<< mol->getNumAtoms() << std::endl;
	}

	Bioisosteres::testCall();

	return 0;

}

