/* 
 * programme to test the implementation in c++ of
 * fingerprints and tanimoto similarity 
 */

// custom headers
#include "Frag.hpp" 
#include "Bioisosteres.hpp"


// rdkit libraries
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

// use RDKit namespace 
using namespace RDKit;

int main()
{
	// insert code here 
	std::string fname = std::string("P39900");
	auto mols = Bioisosteres::getMolsFromFile(fname);
	std::cout<<"There are "<<mols.size()<< " mols in the file." << std::endl;

	// turn two mols into fingerprints 
	
	// declar sparse int vectors for the fingerprints
	
	SparseIntVect<boost::uint32_t> *fp1, *fp2;
	fp1 = MorganFingerprints::getFingerprint(*mols[0],2);
	fp2 = MorganFingerprints::getFingerprint(*mols[1],2);

	// test similarity 
	std::cout<< TanimotoSimilarity(*fp1, *fp1) << std::endl;
	std::cout<< TanimotoSimilarity(*fp2, *fp1) << std::endl;
	
	// fragment a molecule and test the similarity
	auto frags = Bioisosteres::fragment_mol(*mols[0]);

	frags[0].setFingerprint();
	frags[1].setFingerprint();

	std::cout<< TanimotoSimilarity(*frags[0].fp, *frags[1].fp) << std::endl;
	/*	
	for(const auto& fr : frags)
	{
		fr.setFingerprint();
	}
	*/
	delete fp1;
	delete fp2;
	return 0;
}
