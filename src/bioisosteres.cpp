// standard includes
#include <stdlib.h> // exit getenv

// rdkit includes
#include <RDGeneral/Invariant.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/FileParsers/MolSupplier.h> // for obtaining the molecules from sdf
#include <GraphMol/FileParsers/MolWriters.h> // for writing the molecules to sdf
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/ChemTransforms/MolFragmenter.h>
#include <GraphMol/RDKitQueries.h>
#include <DataStructs/SparseIntVect.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>


using namespace RDKit;

typedef boost::shared_ptr<SparseIntVect <boost::uint32_t> > SIV_SPTR;
typedef std::vector<ROMOL_SPTR> MOL_FRAGS;
// typedef boost::shared_ptr<ROMol> ROMOL_SPTR; // defined in chemtransforms

// TODO :: CREATE A STRUCT FOR THE FRAGMENTS THAT WOULD BE ABLE TO HOLD A FRAG A FP AND A GAUSSIAN

void ReadMols(std::string fname, std::vector<ROMOL_SPTR> &mols)
{
// --------------------------------------------------------------
// 	READ MOLECULES
// --------------------------------------------------------------

	char *data = getenv("DATA");
	if (data==NULL)
	{
		std::cerr << "cannot find data environment" << std::endl;
		exit(EXIT_FAILURE);
	}
	

	std::string fpath = std::string(data) + "validation_overlays/" +fname + ".sdf";
	SDMolSupplier suppl(fpath,true);	// sanitize mols and keep H atoms
	while(!suppl.atEnd())
	{
		ROMol *nmol = suppl.next();
		if(!nmol) continue; 			// NULL is 0 in c++
		ROMOL_SPTR p_nmol(nmol);
		mols.push_back(p_nmol);
	}
}


void FragmentMols(const std::vector<ROMOL_SPTR> &mols, std::vector<MOL_FRAGS> &fragments)
{
// --------------------------------------------------------------
// 	CONSTRUCT FRAGMENTS
// --------------------------------------------------------------

	BOOST_FOREACH(ROMOL_SPTR p_mol, mols)
	{
		ROMol* frag = MolFragmenter::fragmentOnBRICSBonds(*p_mol);
		MOL_FRAGS mol = MolOps::getMolFrags(*frag);
		fragments.push_back(mol);
	}
}


void BuildFps(const std::vector<ROMOL_SPTR> &mols, std::vector<SIV_SPTR> &fingerprints)
{
// --------------------------------------------------------------
// 	CONSTRUCT FINGERPRINTS 
// --------------------------------------------------------------

	BOOST_FOREACH(ROMOL_SPTR p_mol, mols)
	{
		SparseIntVect<boost::uint32_t> *fp;
		fp =  MorganFingerprints::getFingerprint(*p_mol,2);
		SIV_SPTR p_fp(fp);
		fingerprints.push_back(p_fp);
	}
}

int main()
{
	std::string file_name = "P39900";
	std::vector<ROMOL_SPTR> mols;
	std::vector<MOL_FRAGS> frags;
	std::vector <SIV_SPTR> mol_fps;
	std::vector<std::vector <SIV_SPTR> > frag_fps;
	std::vector<std::vector <ROMOL_SPTR> > total_sections;

	// get molecules
	std::cout<<"Getting mols... " << std::endl;
	ReadMols(file_name,mols);

	// fragment mols to obtain vectors of fragmented mols
	FragmentMols(mols,frags);

	// get fingerprints of each fragment
	for(unsigned int i = 0; i < frags.size(); ++i)
	{
		std::vector<SIV_SPTR> f_fps;
		BuildFps(frags[i],f_fps);
		frag_fps.push_back(f_fps);
	}

	// print frags and fingerprints to ensure they are the same indices	
	std::cout<<"print frags and fingerprints to ensure they are the same indices" << std::endl;	
	for(unsigned int i = 0; i < frags.size(); ++i)
	{
		std::cout<<"no frags: " << frags[i].size() << "\tno. fps: " << frag_fps[i].size() << std::endl;
	}
	
	int pair_count = 0;

	for(unsigned int ref = 0; ref < frags.size()-1; ++ref)
	{
		MOL_FRAGS ref_mol = frags[ref];

		for(unsigned int query = ref+1; query < frags.size(); ++query)
		{
			MOL_FRAGS query_mol = frags[query];
			// for each fragment in a reference molecule
			// produce the section score
			// get vector of coordinates

			for(unsigned int r_fr = 0; r_fr < ref_mol.size(); ++r_fr)
			{
				ROMOL_SPTR ref_frag = ref_mol[r_fr];

				// create a vector of group pairs
				std::vector<ROMOL_SPTR> section_group;
				section_group.push_back(ref_frag);

				// get the conformer of the ref_frag
				Conformer ref_conf = ref_frag->getConformer();
				auto ref_pos = ref_conf.getPositions();

				
				for(unsigned int q_fr = 0; q_fr < query_mol.size(); ++q_fr)
					{
						// get the conformer of the ref_frag
						ROMOL_SPTR q_frag = query_mol[q_fr];
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
						
						// Calculate Tanimoto Similarity
						double sim = TanimotoSimilarity(*frag_fps[ref][r_fr],*frag_fps[query][q_fr]);
						if (section_score > 0.7 && sim != 1)
						{
							// debugging print out
							//std::cout << "kenn: " << section_score << "\tsim: " << sim << std::endl;
							section_group.push_back(q_frag);
							pair_count++;
						}
					
					}
			
					
				if(section_group.size() > 1)
				{
					std::cout << "number of pairs in section: " << section_group.size() << std::endl;
					total_sections.push_back(section_group);
					/*
					char* home = getenv("HOME");
					string fname = std::string(home) + "fragments/test_output/cpp/pair_" + to_string(count)+".sdf";
					SDWriter *writer = new SDWriter(fname);
					for(auto const& match : section_group)
					{
						writer->write(*match);
					}
					writer->flush();
					writer->close();
					delete writer;
					*/
				}
			}


		}
	}


	std::cout<<"pair count: " << pair_count << std::endl;
	std::cout<<"group count: " << total_sections.size() << std::endl;
		

	return 0;
}
