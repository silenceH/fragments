#ifndef FRAGMENTS
#define FRAGMENTS

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h> // for obtaining the molecules from sdf
#include <GraphMol/FileParsers/MolWriters.h> // for writing the molecules to sdf
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>

using namespace RDKit;
using namespace std;

typedef boost::shared_ptr<ROMol> sp_fragments;

class Fragments
{
	public:
		// ROMol* parent; 	//necessary?? 
		std:: string smiles;
		void generate_fragments(ROMol& mol);
		sp_fragments get_mol_frags();

	
	private:
		ROMol* frag_set;
		sp_fragments frags;


};

#endif
