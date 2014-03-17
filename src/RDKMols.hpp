#ifndef RDKMOLS_DEF

#define RDKMOLS_DEF

/*
 * KENNEWELL_DEV_11.CPP
 * Implementation of Kennewell in c++ using the C++11 standard
 * Matthew Seddon
 * January 2014
 * compile with:
 * g++ -std=c++11 -o frags_11 frags_11.cpp -I$RDBASE/Code -I$RDBASE/EXtern 
 * -L$RDBASE/lib -lFileParsers -lGraphMol -lRDGeneral -lChemTransforms
 */
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h> // exit getenv

#include <GraphMol/ROMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h> // for obtaining the molecules from sdf
#include <GraphMol/FileParsers/MolWriters.h> // for writing the molecules to sdf
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/RDKitQueries.h>

using namespace RDKit;
typedef boost::shared_ptr<ROMol> sp_fragments;


class RDKMols
{
	public:
		// constructor
		RDKMols();
		RDKMols(const std::string file_name);

		// attributes
		std::string ligand;
		std::vector<ROMol*> mols;

		//method
		std::vector<sp_fragments> fragment_mol(ROMol& mol);

		/* 
		 * Method that takes a file name and returns a std::vector of ROMol pointers
		 */
		//void getMols(std::string file_name);

};
#endif
