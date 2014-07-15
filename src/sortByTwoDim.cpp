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

// local files
#include "bioisosteres.cpp"

using namespace RDKit;

typedef boost::shared_ptr<SparseIntVect <boost::uint32_t> > SIV_SPTR;
typedef std::vector<ROMOL_SPTR> MOL_FRAGS;

struct FragmentPair {
        vector<MOL_FRAGS> frags;
        double sim;
};

std::vector<FragmentPair> generateFragmentPairs (std::vector<FRAG>)
{
        std::vector pairVector;
        for(auto const frag1:){
                for(auto const frag2:){
                        FragmentPair add_pair;

                        // add fragments to  pair
                        add_pair.pair.push_back(frag1);
                        add_pair.pair.push_back(frag2);

                        // get similarity
                        add_pair.sim = get_similarity(frag1,frag2);

                        // add pair to vector
                        pairVector.push_back(add_pair);
                }
        }
        return pairVector;
}

                
bool compare_by_similarity (fragment_pair pair1, fragment_pair pair2)
{
        if(pair1.sim >= pair2.sim)
                return true;
}

// sort a vector of pair objects by the similarity

int main()
{
       // get a vector of pairs
        vector<fragment_pair> vector_of_fragment_pairs;
       // sort list
        std::sort(vector_of_fragment_pairs.begin(),vector_of_fragment_pairs.end(),compare_by_similarity);

       for (int i=0;i < vector_of_fragment_pairs.size(); ++i)
               std::cout<<vector_of_fragment_pairs[i].sim;

}
