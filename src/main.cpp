#include <iostream>
#include "Frag.hpp"
#include "FragGroup.hpp"
#include "Bioisosteres.hpp"
#include "RDKMols.hpp"


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
		auto frags = Bioisosteres::fragment_mol(*mol);

	}

	Bioisosteres::testCall();

	// test FragGroup object
	FragGroup* p_group1 = new FragGroup(new_frag4); 

	FragGroup* p_group2 = new FragGroup(*p_group1); 

	std::cout<<"len group1 = " << p_group1->getLength() << std::endl;	
	std::cout<<"len group2 = " << p_group2->getLength() << std::endl;	

	return 0;

}

