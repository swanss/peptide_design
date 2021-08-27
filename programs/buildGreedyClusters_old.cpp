//
//  buildGreedyClusters.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/8/21.
//  Copyright © 2021 Sebastian Swanson. All rights reserved.
//

#include <stdio.h>

#include "mstoptions.h"
#include "msttypes.h"
#include "structure_iter.h"
#include "mstrotlib.h"
#include "mstsystem.h"

string getSeedWindowName(vector<vector<Atom*>>& seed_windows, int window_id, int window_length) {
    Atom* cluster_representative_atom = seed_windows[window_id][0];
    Structure* seed = cluster_representative_atom->getStructure();
    string seed_name = seed->getName();
    int position = cluster_representative_atom->getResidue()->getNum();
    string seed_window_name = seed_name + MstUtils::toString(position) + "-" + MstUtils::toString(position+window_length);
    return seed_window_name;
}

int main(int argc, char* argv[]) {
    MstOptions op;
    op.setTitle("Clusters a set of seeds using a greedy algorithm. For convenience, writes pdb files for the cluster representatives and a few of the cluster members.");
    op.addOption("bin", "Path to a binary file containing seed structures", true);
    op.addOption("cluster_radius", "The RMSD cutoff that is used to identity how many seeds are covered by a representative");
    op.addOption("window_length", "The seed window length that is used when calculating RMSD");
    op.addOption("N_max","If more than this many seeds are present, will be extra greedy");
    op.addOption("coverage","Specifies what fraction of the seeds should be covered in the clusters");
    op.addOption("nClusters","The program will continue until this many clusters have been created");
    op.addOption("pdb","A PDB file describing a peptide-protein complex. If this is provided, will try to map clusters to the peptide");
    op.addOption("peptideChain","The chain ID of the peptide. Required if a complex is provided");
    op.setOptions(argc, argv);
    
    if (op.isGiven("complex") && !op.isGiven("peptideChain")) MstUtils::error("If a complex is provided, a peptide chain ID must also be provided");
    
    string bin_name = op.getString("bin");
    mstreal cluster_RMSD = op.getReal("cluster_radius",1.0);
    int window_length = op.getInt("window_length",5);
    int N_max = op.getInt("N_max",10000);
    mstreal coverage = op.getReal("coverage",1.0);
    int nClusters = op.getReal("nClusters",1000);
    string complex_pdb = op.getString("pdb","");
    string peptideChain = op.getString("peptideChain","");
    
    string seed_chain_id = "0";
    
    bool makeParents = true;
    string outDir = "seed_cluster_members/";
    MstSys::cmkdir(outDir,makeParents);
    
    // Read seeds from file
    StructuresBinaryFile* bin = new StructuresBinaryFile(bin_name);
    long structure_count = bin->structureCount();
    StructureCache* structures = new StructureCache(bin,structure_count);
    structures->preloadFromBinaryFile();
    
    // Determine the subsample rate
    int number_of_seeds_to_write = 1000;
    int number_of_clusters_to_write = 500;
    /*
     This accounts for the fact that there are multiple seed windows per seed. It assumes that on
     average seeds are 15 residues long.
     */
    int adjustment_factor = 15 - window_length + 1; //e.g. expected number of windows
//    mstreal subsample_rate = min(number_of_seeds_to_write/(mstreal(structure_count)*adjustment_factor),1.0);
    mstreal subsample_rate = 1.0;
    
    // Divide seeds into overlapping windows
    vector<vector<Atom*>> seed_windows;
    auto it = structures->begin();
    int count = 1;
    while (it != structures->end()) {
        Structure* extended_fragment = *it;
        Chain* seed_chain = extended_fragment->getChainByID("0");
        // seeds only contain backbone atoms
        vector<Atom*> seed_atoms = seed_chain->getAtoms();
        for (int i = 0; i < seed_chain->residueSize() - window_length + 1; i++) {
            int start = 4*i;
            int end = start+4*window_length;
            vector<Atom*> seed_window = vector<Atom*>(seed_atoms.begin()+start,seed_atoms.begin()+end);
            seed_windows.push_back(seed_window);
        }
        count++;
        it++;
    }
    
    // Perform clustering on seed windows
    bool optimAlignment = false; // the seed structures are fixed in space
    Clusterer cl(optimAlignment);
    cout << "Begin clustering.." << endl;
    vector<vector<int>> clusters = cl.greedyCluster(seed_windows, cluster_RMSD, N_max, coverage, nClusters);
    
    // If a peptide chain is provided, try to map the clusters to the peptide chain
    map<int,pair<int,mstreal>> cluster2peptideWindow; // Maps the cluster ID to the N-terminal residue of the overlapping peptide window
    if (complex_pdb != "") {
        Structure peptideProteinComplex(complex_pdb);
        Structure peptide = *peptideProteinComplex.getChainByID(peptideChain);
        RMSDCalculator rmsdCalc;
        vector<Atom*> peptideAtoms = peptide.getAtoms();
        for (int i = 0; i < clusters.size(); i++) {
            vector<Atom*> cluster_centroid = seed_windows[clusters[i][0]];
            for (int pos = 0; pos < peptide.residueSize(); pos++) {
                vector<Atom*> peptide_window(peptideAtoms.begin()+(pos*4),peptideAtoms.begin()+(pos+window_length)*4);
                mstreal rmsd = rmsdCalc.rmsd(cluster_centroid,peptide_window);
                if (rmsd < cluster_RMSD) {
                    if (cluster2peptideWindow.count(i) == 0) cluster2peptideWindow[i] = pair<int,mstreal>(pos,rmsd);
                    else if (rmsd < cluster2peptideWindow[i].second) cluster2peptideWindow[i] = pair<int,mstreal>(pos,rmsd);
                }
            }
        }
    }
    
    // Write top clusters to PDB files and write out all cluster data
    fstream out;
    MstUtils::openFile(out,"cluster_info.csv",fstream::out);
    out << "cluster_number,cluster_size,cluster_representative,overlapping_peptide_window,rmsd" << endl;
    for (int i = 0; i < clusters.size(); i++) {
        string seed_window_name = getSeedWindowName(seed_windows,clusters[i][0],window_length);
        out << i << ",";
        out << clusters[i].size() << ",";
        out << seed_window_name << ",";
        out << cluster2peptideWindow[i].first << ",";
        out << cluster2peptideWindow[i].second << endl;
    }
    out.close();
    
    for (int i = 0; i < clusters.size(); i++) {
        if (i >= number_of_clusters_to_write) break;
        string cluster_centroid_name = getSeedWindowName(seed_windows,clusters[i][0],window_length);
        Structure cluster_centroid = Structure(seed_windows[clusters[i][0]]);
        string pdb_name = outDir+MstUtils::toString(i)+"_"+MstUtils::toString(0)+"_"+cluster_centroid_name+".pdb";
        cluster_centroid.writePDB(pdb_name);
        for (int j = 1; j < clusters[i].size(); j++) {
            if (MstUtils::randUnit() > subsample_rate) {
                continue;
            }
            string seed_window_name = getSeedWindowName(seed_windows,clusters[i][j],window_length);
            Structure seed_window = Structure(seed_windows[clusters[i][j]]);
            string pdb_name = outDir+MstUtils::toString(i)+"_"+MstUtils::toString(j)+"_"+seed_window_name+".pdb";
            seed_window.writePDB(pdb_name);
        }
    }
}
