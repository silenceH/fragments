// standard libraries
#include <iostream>
#include <vector>

// RDKit headers
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

// custom header files
#include "Bioisosteres.hpp"

using namespace RDKit;

namespace Bioisosteres{
	void testCall()
	{
		std::cout << "called function" << std::endl;
	}

	std::vector<sp_fragments> fragment_mol(ROMol& mol)
	{
		// take a dereferenced mol pointer and fragment
		ROMol* frag = MolFragmenter::fragmentOnBRICSBonds(mol);
		std::vector<sp_fragments> num_frags = MolOps::getMolFrags(*frag);
		return num_frags;
	}

	std::vector<ROMol*> getMolsFromFile(std::string file_name)
	{
		char* data = getenv("DATA");
		if (data==NULL)
		{
			std::cerr << "cannot find data environment" << std::endl;
			exit(EXIT_FAILURE);
		}
		
		// set molecules in file
		std::string fname = std::string(data) + std::string("validation_overlays/") +file_name + std::string(".sdf");
		SDMolSupplier suppl(fname,true);	// sanitize mols and keep H atoms

		std::vector<ROMol *> mols;

		while(!suppl.atEnd())
		{
			ROMol *nmol = suppl.next();
			if(nmol!=0) 			// NULL is 0 in c++
				mols.push_back(nmol);
		}

		return mols;
	}
}
