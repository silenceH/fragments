// sample c++ code taken from rdkit site
// copyright Greg Landrum

#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/RDLog.h>
#include <vector>
#include <algorithm>

using namespace RDKit;

void BuildSimpleMol()
{
	// build the molecule C/C=C\C
	RWMol *mol = new RWMol();

	// add atoms and bonds
	mol->addAtom(new Atom(6)); // atom 0
	mol->addAtom(new Atom(6)); // atom 1
	mol->addAtom(new Atom(6)); // atom 2
	mol->addAtom(new Atom(6)); // atom 3
	mol->addBond(0,1,Bond::SINGLE); // bond 0
	mol->addBond(1,2,Bond::DOUBLE); // bond 1
	mol->addBond(2,3,Bond::SINGLE); // bond 2

	// set up the stereochemistry
	mol->getBondWithIdx(0)->setBondDir(Bond::ENDUPRIGHT);
	mol->getBondWithIdx(2)->setBondDir(Bond::ENDDOWNRIGHT);

	// do chemistry perception
	RDKit::MolOps::sanitizeMol(*mol);

	// get the canonical smiles and include stereochemistry
	std::string smiles;
	smiles=MolToSmiles(*(static_cast<ROMol *>(mol)),true);
	BOOST_LOG(rdInfoLog)<<" sample 1 SMILES: " <<smiles<<std::endl;
}

void WorkWithRingInfo()
{
	// a method for a more complicated job
	ROMol *mol=SmilesToMol("OC1CCC2C1CCCC2");
	// SmilesToMol method is already sanitized so we don't need 
	// to worry about sanitizing the mol this time.
	
	// work with ring information
	RingInfo *ringInfo=mol->getRingInfo();

	TEST_ASSERT (ringInfo->numAtomRings()==2); // returns bool(?)
	
	// ask how many rings an atom is in
	TEST_ASSERT (ringInfo->numAtomRings(0)=0); 
	TEST_ASSERT (ringInfo->numAtomRings(1)==1);
	TEST_ASSERT (ringInfo->numAtomRings(4)==2);


	// the same can be done with bonds
	TEST_ASSERT (ringInfo->numBondRings(0)==0);
	TEST_ASSERT (ringInfo->numBondRings(1)==1);
	
}
int main(int argc, char *argv[])
{
	RDLog::InitLogs();
	BuildSimpleMol();
}




