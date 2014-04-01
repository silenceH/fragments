#ifndef FRAG_DEF
#define FRAG_DEF
//Frag object file

// standard headers
#include <string>

// rdkit headers
#include <GraphMol/ROMol.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>

using namespace RDKit;

typedef boost::shared_ptr<ROMol> sp_fragments;

class Frag
{
	public:
		SparseIntVect<boost::uint32_t> * fp;
		sp_fragments frag;
		// constructor
		Frag();
		Frag(sp_fragments fragment);
		Frag(sp_fragments fragment, std::string file_name);
		
		// set methods
		void setSmiles(std::string finger);
		void setFile(std::string file_name);
		void setFingerprint();

		// get methods
		std::string getSmiles() const;
		std::string getFile() const;

	private:
		// attributes
		std::string smiles;
		std::string file;
};
#endif
