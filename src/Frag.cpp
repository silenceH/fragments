#include "Frag.hpp"

// rdkit headers
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

// constructors
Frag::Frag()
{
	file = "empty";
}

Frag::Frag(sp_fragments fragment)
{
	frag = fragment;
}


Frag::Frag(sp_fragments fragment, std::string file_name)
{
	frag = fragment;
	file = file_name;
}


// set methods
void Frag::setFile(std::string file_name) 
{
	file = file_name;
}

void Frag::setSmiles(std::string s) 
{
	smiles = s;
}

void Frag::setFingerprint()
{
	fp = MorganFingerprints::getFingerprint(*frag,2);
}

// get methods
std::string Frag::getSmiles() const
{
	return smiles;
}

std::string Frag::getFile() const
{
	return file;
}
