#ifndef FRAG_DEF
#define FRAG_DEF
//Frag object file

// standard headers
#include <string>

// rdkit headers
#include <GraphMol/ROMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>

using namespace RDKit;

typedef boost::shared_ptr<ROMol> sp_fragments;

class Frag
{
	public:
		// constructor
		Frag();
		Frag(sp_fragments fragment);
		Frag(sp_fragments fragment, std::string file_name);
		
		// set methods
		void setSmiles(std::string finger);
		void setFile(std::string file_name);

		// get methods
		std::string getSmiles();
		std::string getFile();

	private:
		// attributes
		sp_fragments frag;
		std::string smiles;
		std::string file;
};
#endif
