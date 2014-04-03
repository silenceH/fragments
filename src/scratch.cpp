/*
 * method to return a vector of fragment groups
 */

std::vector<FragGroup> getBioisostericGroups(std::string file_name)
{
	// obtain molecules
	std::vector<ROMol *> mols = Bioisosteres::getMolsFromFile(file_name);
	int num_pairs = 0;
	//std::vector<FragGroup*> section_groups;
	for(vector<ROMol*>::iterator i = mols.begin(); i!=mols.end()-1;++i)
	{
		ROMol *mol = *i;
		std::vector<Frag> ref_fragments = Bioisosteres::fragment_mol(*mol); 	// as mol is a smart pointer we use *mol

		// for the remaining molecules, fragment and score
		for(vector<ROMol*>::iterator j = i+1; j!=mols.end();++j)
		{
			ROMol* query = *j;
			std::vector<Frag> query_fragments = Bioisosteres::fragment_mol(*query);

			// for each fragment in a reference molecule
			// produce the section score
			// get vector of coordinates
			for(const auto& ref_frag : ref_fragments)
			{
				// TODO DO ALL REF_FRAGS HAVE TO BE POINTERS??
				// create a vector of matched pairs
				FragGroup section_match = FragGroup(ref_frag);

				Conformer ref_conf = (*ref_frag.frag)->getConformer();
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
					string fname = "/usr/users/people/matts/fragments/test_output/cpp/pair_" + to_string(num_pairs)+".sdf";
					SDWriter *writer = new SDWriter(fname);
					for(auto const& match : section_match)
					{
						writer->write(*match);
					}
					writer->flush();
					writer->close();
					delete writer;
					num_pairs ++;
				}
			} 
		}
	}
	cout << "number of pairs: " << num_pairs << endl;
	for(const auto& mol : mols)
	{
		delete mol;
	}

}

// NOTE:: we have stuck with the blow code. see frags_11 for other version.
// for(vector<ROMol*>::iterator i = mols.begin(); i!=mols.end()-1;++i)
// for(vector<ROMol*>::iterator j = i+1; j!=mols.end();++j)
