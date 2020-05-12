//
//  scoreSeeds.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 9/12/19.
//

//TO DO

#include "coverage.h"

#include "msttypes.h"
#include "mstoptions.h"
#include "dtermen.h"

#include <regex>

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Takes a set of seeds and the structure of the peptide-protein interaction that they cover and uses the dTERMen framework to generate amino acid probability distributions.");
    op.addOption("structureList","A file containing the paths to all the required structures. The first line is peptide-protein complex and the following lines are the seeds",true);
    op.addOption("peptide","The selection string for the peptide",true);
    op.addOption("c", "dTERMen configuration file.", true);
    op.setOptions(argc, argv);
    
//    ////Import the structures
//    //Import the structure list
//    vector<string> structure_names;
//    structure_names = MstUtils::fileToArray(op.getString("structureList"));
//
//    //Import the complex
//    Structure C(structure_names.front());
//
//    //Import the seeds
//    vector<Structure> seed_structures;
//    for (int i = 1; i < structure_names.size(); i++) {
//        Structure seed(structure_names[i]);
//        seed_structures.push_back(seed);
//    }
//
//    ////Identify the peptide residues that should be scored
//    // Read the structure name to identify the peptide chain name
//    /*
//     PixelDB structures (and in general, any structure used in this benchmark) should be named in the
//     following syntax PDBID_proteinChain_peptideChain, e.g., 5UUL_A_B.
//     */
//    string complex_name(C.getName());
//    regex exp("\S{4}_(\D)_(\D)");
//    smatch match_strings;
//    regex_match(complex_name,match_strings,exp);
//
//    string peptide_chain_name = match_strings[1];
//    cout << "Based on the complex name, the peptide chain is " << MstUtils::toString(peptide_chain_name) << endl;
//
//
//    //Identify the index of the central residue of each seed
//    /*
//     When stored, these seed structures should have been renumbered according to the residue indices
//     of the complex that they cover. These are indices in the overall structure and can be used to
//     identify the corresponding residues in the peptide-protein complex.
//     */
//    vector<int> peptide_res_idx;
//    for (int i = 0; i < seed_structures.size(); i++) {
//        Structure& seed =  seed_structures[i];
//        //NOTE: this could be slow - if it turns out to matter, I will find a different way to get the peptide chain
//        Chain* seed_peptide_chain = seed.getChainByID(peptide_chain_name);
//        //assumes that the first residue in the chain has the lowest residue number
//        int central_res_idx = (*seed_peptide_chain)[0].getNum() + (seed_peptide_chain->residueSize())/2;
//        if (find(peptide_res_idx.begin(), peptide_res_idx.end(), central_res_idx) == peptide_res_idx.end()) {
//            peptide_res_idx.push_back(central_res_idx);
//        }
//    }
//
//    ////Score the residues using the dTERMen statistical framework
//    // Intialize dTERMen-derived class
//    // NOTE: if there are issues, consider extracting the protein first
//    seedScoring SS(op.getString("c"));
//
//    map<int,vector<mstreal>> peptide_res_dist_map;
//    //Score the selected residues in the true peptide-protein complex
//    for (int i = 0; i < peptide_res_idx.size(); i++) {
//        int peptide_res_id = peptide_res_idx[i];
//        vector<mstreal> peptide_res_dist;
//
//        //Build new structure
//        /* So that the probability distributions estimated for the peptide can be comparable to the*/
//
//        //score residue (generate probability distribution)
//        peptide_res_dist = SS.seedDist(&(C.getResidue(peptide_res_id)));
//
//        //add the probability distribution to the map
//        peptide_res_dist_map.emplace(peptide_res_id, peptide_res_dist);
//    }
//
//    //Score the central residues of each seed
//    for (int i = 0; i < seed_structures.size(); i++) {
//
//    }
//
//    //Compute the dot product of the distributions and report results
  
};
