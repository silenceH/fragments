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

struct Frag {
	ROMOL_SPTR frag;
	SIV_SPTR fp;
};

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
	//SDMolSupplier suppl(fpath,true);	// sanitize mols and keep H atoms
	SDMolSupplier suppl(fpath);	// sanitize mols and keep H atoms false
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

void getBioiosteres(std::string file_name, std::vector<std::vector <ROMOL_SPTR> > &final_group, bool debug=false,bool debug_out=false)
{
	std::vector<ROMOL_SPTR> mols;
	std::vector<MOL_FRAGS> frags;
	std::vector <SIV_SPTR> mol_fps;
	std::vector<std::vector <SIV_SPTR> > frag_fps;
	std::vector<std::vector <ROMOL_SPTR> > total_sections;

	// get molecules
	std::cout<<"Getting mols... " << std::endl;
	std::cout<<"file:  " << file_name << std::endl;
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

	// print frags and fingerprints to ensure they are the same indices for debugging
	/*std::cout<<"print frags and fingerprints to ensure they are the same indices" << std::endl;	
	for(unsigned int i = 0; i < frags.size(); ++i)
	{
		std::cout<<"no frags: " << frags[i].size() << "\tno. fps: " << frag_fps[i].size() << std::endl;
	}
	*/
	int pair_count = 0;

	for(unsigned int ref = 0; ref < frags.size()-1; ++ref)
	{
		MOL_FRAGS ref_mol = frags[ref];
		if (debug)
			std::cout << "ref ligand: " << ref << std::endl;

		for(unsigned int query = ref; query < frags.size(); ++query)
		{
			MOL_FRAGS query_mol = frags[query];
			if (debug)
				std::cout << "query ligand: " << query << std::endl;

			// for each fragment in a reference molecule
			// produce the section score
			// get vector of coordinates

			for(unsigned int r_fr = 0; r_fr < ref_mol.size(); ++r_fr)
			{
				if (debug)
					std::cout << "i am reference fragment: " << r_fr << std::endl;

				ROMOL_SPTR ref_frag = ref_mol[r_fr];

				// create a vector of group pairs
				std::vector<ROMOL_SPTR> section_group;
				section_group.push_back(ref_frag);

				// get the conformer of the ref_frag
				Conformer ref_conf = ref_frag->getConformer();
				auto ref_pos = ref_conf.getPositions();

				for(unsigned int q_fr = 0; q_fr < query_mol.size(); ++q_fr)
					{
						if (debug)
							std::cout << "\ti am query fragment: " << q_fr << std::endl;

						// get the conformer of the ref_frag
						ROMOL_SPTR q_frag = query_mol[q_fr];
						Conformer q_conf = q_frag->getConformer();
						auto q_pos = q_conf.getPositions();

						// calculate distance
						double section_score = 0;
						for(const auto& r_atom : ref_pos)
						{
							double atom_score = 0;
							for(const auto& q_atom : q_pos)
							{
								double x = r_atom.x - q_atom.x;
								double y = r_atom.y - q_atom.y;
								double z = r_atom.z - q_atom.z;
								double dist = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
								// calculate gaussian overlap 
								atom_score += exp(-pow(dist,2));
							}
							section_score += atom_score;
						}

						// calculate average overlap 
						double total_atoms = (ref_frag->getNumAtoms() + q_frag->getNumAtoms());
						double av_score = section_score * (2/total_atoms);

						if (debug)
							std::cout << "\taverage score: " << std::setprecision(16) << av_score<<std::endl;
						
						// Calculate Tanimoto Similarity
						double sim = TanimotoSimilarity(*frag_fps[ref][r_fr],*frag_fps[query][q_fr]);
						if (av_score > 0.7 && sim != 1)
						{
							section_group.push_back(q_frag);
							pair_count++;
							if(debug_out){
								std::cout<<"pair: " << pair_count << std::endl;
								std::cout<<"Ref: " << ref << "\tfrag: "<<r_fr<<std::endl;
								std::cout<<"Query: " << query << "\tfrag: "<<q_fr<<std::endl;
							}
						}
					
					}
					
				if(section_group.size() > 1)
				{
					//std::cout << "number of pairs in section: " << section_group.size() << std::endl;
					total_sections.push_back(section_group);

					// merge groups that have molecules in common. 

					bool in_final_group = false;
					for(auto& final_group_section : final_group)
					{
						bool in_group = false;
						BOOST_FOREACH(ROMOL_SPTR p_mol, final_group_section)
						{
							if(section_group[0] == p_mol)
								in_group = true;
						}

						if(in_group)
						{
							in_final_group == true;
							for(int i=1; i<section_group.size(); ++i)
								final_group_section.push_back(section_group[i]);
						}

					}

					if(!in_final_group)
						final_group.push_back(section_group);
								


					char* home = getenv("HOME");
					std::string fname = std::string(home) + "/fragments/test_output/cpp/pair_" + std::to_string(pair_count)+".sdf";
					SDWriter *writer = new SDWriter(fname);
					for(auto const& match : section_group)
					{
						writer->write(*match);
					}
					writer->flush();
					writer->close();
					delete writer;
				}
			}


		}
	}


	std::cout<<"pair count: " << pair_count << std::endl;
	std::cout<<"group count: " << total_sections.size() << std::endl;
	std::cout<<"final group count: " << final_group.size() << std::endl;
	std::cout<<"\n\n";

}

int main()
{
	std::string file_name_1 = "P39900";
	//std::string file_name_2 = "P56817";
	//std::string file_name_3 = "P35557";
	//std::string file_name_4 = "Q92731";
	//std::string file_name_5 = "P25440";
	
	std::vector<std::vector <ROMOL_SPTR> > final_group;

	getBioiosteres(file_name_1,final_group);
	//getBioiosteres(file_name_2,final_group);
	//getBioiosteres(file_name_3,final_group);
	//getBioiosteres(file_name_4,final_group);
	//getBioiosteres(file_name_5,final_group);
}
