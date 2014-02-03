/*
 * KENNEWELL_DEV_11.CPP
 * Implementation of Kennewell in c++ using the C++11 standard
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
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>
#include <RDGeneral/RDLog.h>


using namespace RDKit;
using namespace std;


auto fragment_mol(ROMol& mol)
{
	// focus on BRICS fragmentation with no overlapping
	// ROMol *pattern=SmartsToMol("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]");
	// ROMol* frags = MolFragmenter::fragmentOnBRICSBonds(mol);
	// Introduce substructure match 
	// MatchVectType matchV;
	// vector<MatchVectType> matches;
	// unsigned int n = SubstructMatch(mol,*pattern,matches);
	ROMol* frag = MolFragmenter::fragmentOnBRICSBonds(mol);
	auto num_frags = MolOps::getMolFrags(*frag);
	return num_frags;
}

int[] get_coords(ROMol& mol)
{
	// takes in a mol object and returns the coordinates of that mol
	int[] coords;
	// get conformer object and return coords
	return coords;
}
vector<RMol*> getMols(string file_name)
{
	string fname = string("../data/validation_overlays/") +file_name + string(".sdf");
	SDMolSupplier suppl(fname,true,false);	// sanitize mols and keep H atoms
	vector<ROMol*> mols;
	while(!suppl.atEnd())
	{
		ROMol *nmol = suppl.next();
		if(nmol!=0) 			// NULL is 0 in c++
			mols.push_back(nmol);
	}
	return mols;
}

int main()
{

	auto mols = getMols("P39900");
	cout<<"list size is: " << mols.size() << endl;
	for(int i=0;i<mols.size();i++)
	{
		cout << mols[i]->getNumAtoms() << endl;
	}

	cout << "number of atoms in each fragment: "<< endl;
	auto new_frags = fragment_mol(*mols[0]);
	for(const auto& f : new_frags)
	{
		cout << f->getNumAtoms() << endl;
	}
	

	return 0;

}

// TODO: Fragment Mols
// TODO: Score Mols
// TODO: I/O
