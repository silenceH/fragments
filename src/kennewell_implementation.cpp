/*
 * KENNEWELL_IMPLEMENTATION.CPP
 * Implementation of Kennewell in c++ 
 * Matthew Seddon
 * January 2014
 * compile with:
 * g++ -o sample.exe sample.cpp -I$RDBASE/Code -I$RDBASE/Extern
 * -L$RDBASE/lib -lFileParsers -lDepictor -lChemTransforms
 * -lSubstructMatch -lGraphMol -lDataStructs -lRDGeometryLib -lRDGeneral
 */
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/Invariant.h>
#include <GraphMol/FileParsers/MolSupplier.h> // for obtaining the molecules from sdf
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <RDGeneral/RDLog.h>


using namespace RDKit;
using namespace std;


vector<boost::shared_ptr<ROMol> > fragment_mol(ROMol& mol)
{
	// ROMol *pattern=SmartsToMol("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]");
	// ROMol* frags = MolFragmenter::fragmentOnBRICSBonds(mol);
	ROMol* frag = MolFragmenter::fragmentOnBRICSBonds(mol);
	vector<boost::shared_ptr<ROMol> > num_frags = MolOps::getMolFrags(*frag);
	return num_frags;
}

int main()
{
	string fname = "../data/validation_overlays/P39900.sdf";
	SDMolSupplier suppl(fname,true,false);	// keep H atoms
	vector<ROMol*> mols;
	while(!suppl.atEnd())
	{
		ROMol *nmol = suppl.next();
		if(nmol!=0) 			// NULL is 0 in c++
			mols.push_back(nmol);
	}

	cout<<"list size is: " << mols.size() << endl;
	for(int i=0;i<mols.size();i++)
	{
		cout << mols[i]->getNumAtoms() << endl;
	}

	cout << "number of atoms in each fragment: "<< endl;
	vector<boost::shared_ptr<ROMol> > new_frags = fragment_mol(*mols[0]);
	for(int i=0;i<new_frags.size();i++)
	{
		cout << new_frags[i]->getNumAtoms() << endl;
	}
	
	return 0;

}

