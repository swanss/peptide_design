//
//  fuser.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 1/15/19.
//

#include <stdio.h>
#include "fusehelper.h"

void FuseHelper::writeFuseResultsFromFile(FuseCandidateFile file, string outDir, string pdbPrefix) {
    MstUtils::assert(outDir.length() > 0 && outDir.back() == '/', "output directory must end with /");
    
    // Output CSV file
    SeedListFile summaryFile(outDir + "fuse_results.csv");
    summaryFile.setMetadataFields({ "total_score", "bond_score", "angle_score", "dihedral_score", "rmsd_score", "total_rmsd_score", "source_1", "source_chain_1", "source_2", "source_chain_2" });
    
    vector<FuseCandidate> candidateBatch;
    int fuseIndex = 0;
    while ((candidateBatch = file.read(1)).size() > 0) {
        FuseCandidate candidate = candidateBatch[0];
        candidate.loadStructures(pdbPrefix);

        string pdbName = "fuse_result_" + to_string(fuseIndex++) + ".pdb";
        tuple<Structure, string, fusionOutput> result = performFuse(candidate);
        get<0>(result).writePDB(outDir + pdbName);
        fusionOutput output = get<2>(result);

        summaryFile.write(pdbName, get<1>(result), {
            // Scores to write along with the seed info
            to_string(output.getScore()),
            to_string(output.getBondScore()),
            to_string(output.getAngleScore()),
            to_string(output.getDihedralScore()),
            to_string(output.getRMSDScore()),
            to_string(output.getTotRMSDScore()),
            candidate.file1,
            candidate.chain1,
            candidate.file2,
            candidate.chain2
        });
        candidate.freeStructures();
    }
    cout << "Fused " << fuseIndex << " candidates" << endl;
}

tuple<Structure, string, fusionOutput> FuseHelper::performFuse(FuseCandidate candidate) {
    cout << candidate.structure1->getName() << ", " << candidate.structure2->getName() << endl;
    
    vector<vector<Residue *>> residues(candidate.structure1->residueSize() + candidate.structure2->residueSize() - candidate.overlapSize);
    // Add all residues not in the overlapping chains
    int residueIdx = 0;
    for (int i = 0; i < candidate.structure1->residueSize(); i++) {
        Residue &res = candidate.structure1->getResidue(i);
        if (res.getChain()->getID() != candidate.chain1) {
            residues[residueIdx++].push_back(&res);
        }
    }
    for (int i = 0; i < candidate.structure2->residueSize(); i++) {
        Residue &res = candidate.structure2->getResidue(i);
        if (res.getChain()->getID() != candidate.chain2) {
            residues[residueIdx++].push_back(&res);
        }
    }
    
    // Add seed chain residues - chain1 will be at N-terminus
    Chain *chain1 = candidate.structure1->getChainByID(candidate.chain1);
    vector<Residue *> residues1 = chain1->getResidues();
    int seedChainStartIdx = residueIdx;
    for (int i = 0; i < candidate.overlapPosition1 + candidate.overlapSize; i++) {
        residues[residueIdx++].push_back(residues1[i]);
    }
    residueIdx -= candidate.overlapSize;
    Chain *chain2 = candidate.structure2->getChainByID(candidate.chain2);
    vector<Residue *> residues2 = chain2->getResidues();
    for (int i = candidate.overlapPosition2; i < residues2.size(); i++) {
        residues[residueIdx++].push_back(residues2[i]);
    }
    
    if (residueIdx != residues.size())
        cerr << "Missed the mark: expected " << residues.size() << " residues in topology, but got " << residueIdx << endl;
    
    fusionTopology topology(residues);
    for (int i = 0; i < seedChainStartIdx; i++) topology.addFixedPosition(i); // Fix the initial residues
    
    // Fuse!
    fusionOutput output;
    Fuser myFuser;
    Structure fusedStruct = myFuser.fuse(topology, output);
    cout << "Output: " << output.getScore() << endl;
    
    // Figure out which chains aren't fixed
    string seedChains = "";
    for (int i = 0; i < topology.numChains(); i++) {
        MstUtils::assert(topology.getChainLengths()[i] == fusedStruct.getChain(i).residueSize(), "chain mismatch");
        if (topology.numFixedInChain(i) == 0) {
            seedChains += fusedStruct.getChain(i).getID();
        }
    }
    
    return make_tuple(fusedStruct, seedChains, output);
}
