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
class Fragment
{
	// insert code
	public:
		sp_fragments frag;
		string smiles;
		//ExplicitBitVect *fp;

};

/*
 * Method to fragment a mol into BRICS fragments
 */
vector<sp_fragments> fragment_mol(ROMol& mol)
{
	// take a dereferenced mol pointer and fragment
	ROMol* frag = MolFragmenter::fragmentOnBRICSBonds(mol);
	vector<sp_fragments> num_frags = MolOps::getMolFrags(*frag);
	return num_frags;
}


/* 
 * Method that takes a file name and returns a vector of ROMol pointers
 */
vector<ROMol*> getMols(string file_name)
{
	char* data = getenv("DATA");
	if (data==NULL)
	{
		cerr << "cannot find data environment" << endl;
		exit(EXIT_FAILURE);
	}
	

	string fname = string(data) + string("validation_overlays/") +file_name + string(".sdf");
	SDMolSupplier suppl(fname,true);	// sanitize mols and keep H atoms
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
	// get the molecules as a vector of pointers
	// for each mol in mols, make the mol the reference and the 
	// rest the queries
	
	string var = "P39900";
	auto mols = getMols(var);
	int count = 0;
	for(vector<ROMol*>::iterator i = mols.begin(); i!=mols.end();++i)
	{
		ROMol *mol = *i;
		auto ref_fragments = fragment_mol(*mol); 	// as mol is a smart pointer we use *mol

		// for the remaining molecules, fragment and score
		for(vector<ROMol*>::iterator j = i; j!=mols.end();++j)
		{
			if(*j != mol)
			{
				ROMol* query = *j;
				auto query_fragments = fragment_mol(*query);

				// for each fragment in a reference molecule
				// produce the section score
				// get vector of coordinates
				for(const auto& ref_frag : ref_fragments)
				{
					// create a vector of matched pairs
					vector<sp_fragments> section_match;
					section_match.push_back(ref_frag);

					Conformer ref_conf = ref_frag->getConformer();
					auto ref_pos = ref_conf.getPositions();
					for(const auto& q_frag : query_fragments)
					{
						Conformer q_conf = q_frag->getConformer();
						auto q_pos = q_conf.getPositions();

						// calculate distance
						long double section_dist = 0;
						for(const auto& r_atom : ref_pos)
						{
							for(const auto& q_atom : q_pos)
							{
								double x = r_atom.x - q_atom.x;
								double y = r_atom.y - q_atom.y;
								double z = r_atom.z - q_atom.z;
								long double dist = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
								// calculate gaussian overlap 
								section_dist += exp(-pow(dist,2));
							}
						}
						// calculate average overlap 
						double total_atoms = (ref_frag->getNumAtoms() + q_frag->getNumAtoms());
						long double section_score = section_dist * (2/total_atoms);
						// test:: print average overlap if successfull otherwise print "NOT A PAIR" 
						if (section_score > 0.7)
						{
							section_match.push_back(q_frag);
						}
					
					}
					if(section_match.size() > 1)
					{
						//cout << "number of pairs in section: " << section_match.size() << endl;
						string fname = "/usr/users/people/matts/fragments/test_output/cpp/pair_" + to_string(count)+".sdf";
						SDWriter *writer = new SDWriter(fname);
						for(auto const& match : section_match)
						{
							writer->write(*match);
						}
						writer->flush();
						writer->close();
						delete writer;
						count ++;
					}
				}
			} 
		}
	}
	cout << "number of pairs: " << count << endl;
	for(const auto& mol : mols)
	{
		delete mol;
	}

}

// if we assume that all mols are fragmented equally then we only need to
// overlay them once?? 
// for(vector<ROMol*>::iterator i = mols.begin(); i!=mols.end()-1;++i)
// for(vector<ROMol*>::iterator j = i+1; j!=mols.end();++j)
// TODO:: optimise loops?? 
// TODO:: delete fingerprints
