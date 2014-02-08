/* 
 * Test to see how to fragment mols ans return a vector of mols
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include <RDGeneral/utils.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h> // for obtaining the molecules from sdf
#include <GraphMol/FileParsers/MolWriters.h> // for writing the molecules to sdf
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Substruct/SubstructUtils.h>

using namespace RDKit;
using namespace std;

typedef boost::shared_ptr<ROMol> sp_fragments;

int main()
{
	// get molecule
	string smiles = "CCCOCCC(=O)c1ccccc1";
	RWMol *mol = SmilesToMol(smiles);

	cout <<(mol->getNumAtoms() == 14) << endl;

	// fragment the molecule
	ROMol *nmol = MolFragmenter::fragmentOnBRICSBonds(*mol);

	vector<sp_fragments> fragments = MolOps::getMolFrags(*nmol);
	
	cout <<"number of fragments: " << fragments.size() << endl;

	vector<sp_fragments> candidate_pairs;

	for(vector<sp_fragments>::iterator i = fragments.begin();i!=fragments.end(); ++i)
	{
		int num = (*i)->getNumAtoms();
		cout << num << endl;
		if(num == 3 || num == 6)
		{
			candidate_pairs.push_back(*i);
		}
	}
	
	
	cout<<"number of pairs: "<<candidate_pairs.size()<<endl;
	for(vector<sp_fragments>::iterator i = candidate_pairs.begin();i!=candidate_pairs.end(); ++i)
	{
		int num = (*i)->getNumAtoms();
		cout << num << endl;
	}
	
	delete nmol;
	delete mol;
	return 0;
}
