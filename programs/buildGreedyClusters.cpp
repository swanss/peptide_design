#include <stdio.h>

#include "mstoptions.h"
#include "msttypes.h"
#include "structure_iter.h"
#include "mstrotlib.h"
#include "mstsystem.h"
#include "clusterutils.h"

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
    op.setTitle("Clusters a set of seeds using a greedy algorithm. For convenience, writes pdb files for the cluster representatives and the cluster members.");
    op.addOption("overlaps","The path to an overlaps file OR a directory containing overlaps files",true);
    op.addOption("bin", "Path to a binary file containing seed structures", true);
    op.addOption("window_length", "The seed window length that is used when calculating RMSD w.r.t the peptide");
    op.addOption("coverage","Specifies what fraction of the seeds should be covered in the clusters");
    op.addOption("nClusters","The program will continue until this many clusters have been created");
    op.addOption("pdb","A PDB file describing a peptide-protein complex. If this is provided, will try to map clusters to the peptide");
    op.addOption("peptideChain","The chain ID of the peptide. Required if a complex is provided");
    op.setOptions(argc, argv);
    
    if (op.isGiven("complex") && !op.isGiven("peptideChain")) MstUtils::error("If a complex is provided, a peptide chain ID must also be provided");
    
    string overlaps = op.getString("overlaps");
    string bin_name = op.getString("bin");
    int window_length = op.getInt("window_length",4);
    mstreal coverage = op.getReal("coverage",1.0);
    string complex_pdb = op.getString("pdb","");
    string peptideChainID = op.getString("peptideChain","");
        
    bool makeParents = true;
    string outDir = "seed_cluster_members/";
    MstSys::cmkdir(outDir,makeParents);
    
    cout << "Construct the greedy clusterer..." << endl;
    GreedyClusterer clusterer(bin_name,window_length);
    
    cout << "Loading overlaps..." << endl;
    if (overlaps.size() > 0) {
        // Loading overlaps from explicit overlap paths
        if (MstSys::isDir(overlaps)) {
            int i = 1;
            string path = MstSystemExtension::join(overlaps, "overlaps" + to_string(i++) + ".csv");
            cout << "path: " << path << endl;
            while (MstSystemExtension::fileExists(path)) {
                FuseCandidateFile file(path);
                clusterer.addOverlapInfo(file);
                path = MstSystemExtension::join(overlaps, "overlaps" + to_string(i++) + ".csv");
                cout << "path: " << path << endl;
            }
        } else {
            // All overlaps stored in one file
            FuseCandidateFile file(overlaps);
            clusterer.addOverlapInfo(file);
        }
    } else {
        MstUtils::error("Could not intepret overlaps path");
    }
    cout << "Done loading overlaps" << endl;

    clusterer.performClustering(coverage);
    
    Chain* peptideChain = nullptr;
    Structure* complex = nullptr;
    if (complex_pdb != "") {
        Structure* complex = new Structure(complex_pdb);
        if (peptideChainID != "") {
            peptideChain = complex->getChainByID(peptideChainID);
            if (peptideChain != nullptr) cout << "Selected chain with ID: " << peptideChain->getID() << endl;
        }
    }
    
    cout << "Write cluster info..." << endl;
    clusterer.writeClusterInfo(outDir,peptideChain,false,true);
    
    delete complex;
    
    cout << "Done" << endl;
    return 0;
}
