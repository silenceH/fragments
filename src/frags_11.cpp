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

/*
 * Method to fragment a mol into BRICS fragments
 */
vector<boost::shared_ptr<ROMol> > fragment_mol(ROMol& mol)
{
	// take a dereferenced mol pointer and fragment
	ROMol* frag = MolFragmenter::fragmentOnBRICSBonds(mol);
	auto num_frags = MolOps::getMolFrags(*frag);
	return num_frags;
}


/* 
 * Method that takes a file name and returns a vector of ROMol pointers
 */
vector<ROMol*> getMols(string file_name)
{
	string fname = string("/home/matthew/data/validation_overlays/") +file_name + string(".sdf");
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
	// get the molecules as a vector of pointers
	auto mols = getMols("P39900");
	// for each mol in mols, make the mol the reference and the 
	// rest the queries
	for(vector<ROMol*>::iterator i = mols.begin(); i!=mols.end();++i)
	{
		ROMol* mol = *i;
		cout <<"size of vector: " << mols.size() << endl;


		// for the remaining molecules, fragment and score
		int count = 0;
		for(vector<ROMol*>::iterator j = mols.begin(); j!=mols.end();++j)
		{
			if(*j != mol)
			{
				ROMol* query = *j;
				auto query_fragments = fragment_mol(*query);

				// for each fragment in a reference molecule
				// produce the section score
				// get vector of coordinates
				Conformer ref_conf = mol->getConformer();
				auto ref_pos = ref_conf.getPositions();
				Conformer q_conf = query->getConformer();
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
						section_dist += exp(-pow(dist,2));
					}
				}
				double tot = (mol->getNumAtoms() + query->getNumAtoms());
				long double section_score = section_dist * (2/tot);
				if (section_score > 0.7)
				{
					cout << section_dist << endl;
				} else
				{
					cout << "NOT A PAIR" << endl;
				}

				count++;
			} 
		}
		cout << "number of comparisons: " << count << endl;
	

	}
	auto new_frags = fragment_mol(*mols[0]);

	return 0;

}

// if we assume that all mols are fragmented equally then we only need to
// overlay them once?? 
// for(vector<ROMol*>::iterator i = mols.begin(); i!=mols.end()-1;++i)	optimised ?? 
// for(vector<ROMfor(vector<ROMol*>::iterator j = i+1; j!=mols.end();++j)
