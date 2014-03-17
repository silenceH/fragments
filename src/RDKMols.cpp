#include <vector>
#include <stdlib.h> // exit getenv
#include "RDKMols.hpp"

#include <GraphMol/ROMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/FileParsers/MolSupplier.h> // for obtaining the molecules from sdf
#include <GraphMol/FileParsers/MolWriters.h> // for writing the molecules to sdf
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/RDKitQueries.h>

using namespace RDKit;

/* 
 * Constructor
 */

RDKMols::RDKMols()
{
	ligand = "unspecified";
}

RDKMols::RDKMols(const std::string file_name)
{
	// set file name
	ligand = file_name;

	char* data = getenv("DATA");
	if (data==NULL)
	{
		std::cerr << "cannot find data environment" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// set molecules in file
	std::string fname = std::string(data) + std::string("validation_overlays/") +file_name + std::string(".sdf");
	SDMolSupplier suppl(fname,true);	// sanitize mols and keep H atoms
	while(!suppl.atEnd())
	{
		ROMol *nmol = suppl.next();
		if(nmol!=0) 			// NULL is 0 in c++
			mols.push_back(nmol);
	}
}


