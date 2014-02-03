#include <GraphMol/RDKitBase.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include <RDGeneral/Invariant.h>
#include <GraphMol/FileParsers/MolSupplier.h> // for obtaining the molecules from sdf
#include <GraphMol/FileParsers/FileParsers.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/RDLog.h>


using namespace RDKit;
using namespace std;


int main()
{
	string fname = "Q00511.sdf";
	SDMolSupplier suppl (fname);
	vector<ROMol*> mols;
	while(!suppl.atEnd())
	{
		ROMol *nmol = suppl.next();
		if(nmol!=0) 			// NULL is 0 in c++
			mols.push_back(nmol);
	}

	cout<<"list size is: " << mols.size() << endl;
	for(int i=0;i<mols.size();i++)
	{
		cout << mols[i]->getNumAtoms() << endl;
	}
	return 0;

}

