#include <iostream>
#include "Frag.hpp"
#include "FragGroup.hpp"
#include "Bioisosteres.hpp"
#include "Bioisosteres.hpp"


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
	std::vector<ROMol *> mols = Bioisosteres::getMolsFromFile(std::string("P39900"));

	std::cout << "Number of molecules: " << mols.size() << std::endl;
	

	for(const auto& mol : mols)
	{
		auto frags = Bioisosteres::fragment_mol(*mol);
		FragGroup* p_group_temp = new FragGroup(frags);
		std::cout<<"no of frags: " << p_group_temp->getLength() << std::endl;
	}

	Bioisosteres::testCall();

	// test FragGroup object
	FragGroup* p_group1 = new FragGroup(new_frag4); 


	std::cout<<"len group1 = " << p_group1->getLength() << std::endl;	

	return 0;

}

