//
//  testFusePaths.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 5/29/19.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "findpaths.h"
#include "utilities.h"

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    string graphFile = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/1DUY_overlaps_dot_75/bond_graphs_stringent/whole_graph.txt";
    string seedPath = "/home/ifsdata/scratch/grigoryanlab/swans/term_based_design_method/old_seed_dirs/1DUY/output/seeds/";
    string scoresPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/1DUY_overlaps_dot_75/seed_scores/";
    string targetPath = "/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb";
    string pdbWritePath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/1DUY_overlaps_dot_75/fused_paths_stringent/";
    vector<string> paths = MstUtils::fileToArray("/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/1DUY_overlaps_dot_75/stringent_paths.txt");
    
    SeedGraph graph;
    graph.read(graphFile, seedPath);
    cout << graph.residueSize() << " residues" << endl;
    
    Structure target(targetPath);
    vector<Residue *> targetResidues = target.getResidues();
    
    string fasstDB = "/home/grigoryanlab/library/databases/dTERMen.databases/2019-01-22/dtermen.sim";
    int numSeedFlank = 2; //opts.getInt("numSeedFlank", 2);
    int numTargetFlank = 2; // opts.getInt("numTargetFlank", 2);
    
    // Initialize
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceCompatibilityScorer scorer(&target, rParams, cParams, fasstDB, numTargetFlank, numSeedFlank, 0.4, 0.05, 0.25, 1, 8000, 0.7);
    int fusionIndex = 0;

    for (string line: paths) {
        vector<string> comps = splitString(line, ",");
        vector<vector<Residue *>> residues;
        for (int i = 1; i < comps.size(); i++) {
            if (comps[i].size() == 0)
                continue;
            Residue *seedRes = graph.getResidueFromFile(comps[i]);
            if (seedRes == nullptr)
                continue;
            // Figure out what match this seed residue came from
            int matchNumber = getTargetResidueIndex(seedRes->getStructure()->getName());
            cout << "Match number: " << matchNumber << endl;
            for (int targetIdx = max(matchNumber - 2, 0); targetIdx <= min(matchNumber + 2, target.residueSize() - 1); targetIdx++) {
                residues.push_back({ targetResidues[targetIdx] });
            }
        }
        
        int seedChainStartIdx = residues.size();
        for (int i = 1; i < comps.size(); i++) {
            if (comps[i].size() == 0)
                continue;
            Residue *seedRes = graph.getResidueFromFile(comps[i]);
            if (seedRes == nullptr)
                continue;
            residues.push_back({ seedRes });
        }
        
        // Fuse the path
        fusionTopology topology(residues);
        for (int i = 0; i < seedChainStartIdx; i++) topology.addFixedPosition(i); // Fix the initial residues
        
        fusionOutput output;
        Fuser myFuser;
        Structure fusedStruct = myFuser.fuse(topology, output);
        // Figure out which chains aren't fixed
        Structure fusedSeed;
        for (int i = 0; i < topology.numChains(); i++) {
            MstUtils::assert(topology.getChainLengths()[i] == fusedStruct.getChain(i).residueSize(), "chain mismatch");
            if (topology.numFixedInChain(i) == 0) {
                fusedSeed.appendChain(new Chain(fusedStruct.getChain(i)));
            }
        }
        fusedSeed.writePDB(MstSystemExtension::join(pdbWritePath, "fusion_" + to_string(fusionIndex) + ".pdb"));
        
        unordered_map<Residue *, mstreal> scoreMap = scorer.score(&fusedSeed);
        double totalScore = 0.0;
        for (auto item: scoreMap) {
            totalScore += item.second;
        }
        
        cout << "Output: " << output.getScore() << endl;
        cout << "Original score: " << comps[0] << endl;
        cout << "New score: " << totalScore << endl;
        fusionIndex++;
    }
    
    return 0;
}
