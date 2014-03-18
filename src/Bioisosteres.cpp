#include "Bioisosteres.hpp"
#include <iostream>
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
}
