//
//  findOverlaps.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 3/22/19.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "overlaps.h"
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Finds overlaps between all pairs of the given seed residues.");
    //opts.addOption("files", "Text file containing a list of seed structures", true);
    //opts.addOption("data", "Directory in which seed structures are stored", true);
    opts.addOption("bin", "Path to a binary file containing seed structures", true);
    opts.addOption("out", "Path to CSV file at which to write overlaps", true);
    opts.addOption("overlapSize", "Number of residues that must overlap between two residues (default 2)", false);
    opts.addOption("overlapRMSD", "Maximum RMSD between sets of residues to be considered overlapping (default 1.0)", false);
    opts.addOption("minCosAngle", "Minimum cosine angle between Ca vectors allowed in a potential overlap (default -1.0)", false);
    opts.addOption("batch", "The batch number (from 1 to numBatches)", false);
    opts.addOption("numBatches", "The number of batches", false);
    opts.addOption("batchSize", "The number of structures to use in each batch (default 200,000)", false);
    opts.addOption("bruteForce", "If provided, use a brute-force all-to-all comparison", false);
    opts.addOption("mock", "If provided, do not perform any overlap calculations, just list the batches that would be computed", false);
    opts.setOptions(argc, argv);
    
    //string filesPath = opts.getString("files");
    //string dataPath = opts.getString("data");
    string binaryPath = opts.getString("bin");
    string outPath = opts.getString("out");

    int numResOverlap = opts.getInt("overlapSize", 2);
    float rmsdCutoff = opts.getReal("overlapRMSD", 1.0);
    float minCosAngle = opts.getReal("minCosAngle", -1.0);

    int batchIndex = opts.getInt("batch", 1);
    int numBatches = opts.getInt("numBatches", 1);
    int batchSize = opts.getInt("batchSize", 200000);
    if (batchIndex < 1 || batchIndex > numBatches) {
        cerr << "Batch index must be between 1 and numBatches" << endl;
        return 1;
    }
    cout << "Batch " << batchIndex << " of " << numBatches << endl;
    
    MaxDeviationVerifier verifier(rmsdCutoff);
    FuseCandidateFile outFile(outPath);
    BatchPairStructureIterator structureIter(binaryPath, batchIndex - 1, numBatches, batchSize); // large batch size
    
    while (structureIter.hasNext()) {
        auto batch = structureIter.next();
        
        if (opts.isGiven("mock")) {
            cout << "Batch: " << batch.first[0]->getName() << ", " << batch.second[0]->getName() << endl;
            continue;
        }
        
        if (opts.isGiven("bruteForce")) {
            // All-to-all comparison
            for (int i = 0; i < batch.first.size(); i++) {
                if (i % 100 == 0)
                    cout << "Structure " << i << " of " << batch.first.size() << endl;
                
                Structure *s1 = batch.first[i];
                
                Chain *c1 = s1->getChainByID("0");
                if (!c1)
                    continue;
                
                vector<FuseCandidate> results;
                
                // Iterate over all segments in the chain
                vector<Residue *> residues1 = c1->getResidues();
                for (int resIdx1 = 0; resIdx1 < residues1.size() - numResOverlap + 1; resIdx1++) {
                    vector<Residue *> segment1(residues1.begin() + resIdx1, residues1.begin() + resIdx1 + numResOverlap);
                    
                    // Compare to all structures
                    for (Structure *s2: batch.second) {
                        if (s2 <= s1) continue;
                        
                        Chain *c2 = s2->getChainByID("0");
                        if (!c2)
                            continue;
                        
                        vector<Residue *> residues2 = c2->getResidues();
                        for (int resIdx2 = 0; resIdx2 < residues2.size() - numResOverlap + 1; resIdx2++) {
                            vector<Residue *> segment2(residues2.begin() + resIdx2, residues2.begin() + resIdx2 + numResOverlap);
                            
                            if (verifier.verify(segment1, segment2)) {
                                FuseCandidate fuseCandidate;
                                fuseCandidate.overlapSize = numResOverlap;
                                fuseCandidate.setStructure1(s1, c1->getID());
                                fuseCandidate.overlapPosition1 = resIdx1;
                                fuseCandidate.setStructure2(s2, c2->getID());
                                fuseCandidate.overlapPosition2 = resIdx2;
                                fuseCandidate.rmsd = 0.0;
                                results.push_back(fuseCandidate);
                            }
                        }
                    }
                }
                
                outFile.write(results, "");
            }
        } else {
            mstreal xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9, zlo = 1e9, zhi = -1e9;
            vector<Structure *> &firstStructures = batch.first;
            vector<Structure *> &secondStructures = batch.second;
            
            auto it = firstStructures.begin();
            while (it != secondStructures.end()) {
                Structure *seed = *it;
                if (seed->residueSize() == 0)
                    continue;
                
                mstreal ixlo, ixhi, iylo, iyhi, izlo, izhi;
                ProximitySearch::calculateExtent(*seed, ixlo, iylo, izlo, ixhi, iyhi, izhi);
                xlo = min(xlo, ixlo);
                xhi = max(xhi, ixhi);
                ylo = min(ylo, iylo);
                yhi = max(yhi, iyhi);
                zlo = min(zlo, izlo);
                zhi = max(zhi, izhi);
                
                ++it;
                if (it == firstStructures.end())
                    it = secondStructures.begin();
            }
            cout << "Bounding box: " << xlo << ", " << xhi << ", " << ylo << ", " << yhi << ", " << zlo << ", " << zhi << ", " << endl;
            vector<mstreal> bbox = { xlo, xhi, ylo, yhi, zlo, zhi };
            
            CAResidueHasher<> hasher(bbox, 0.3);
            OverlapFinder<CAResidueHasher<>> overlapFinder(hasher, rmsdCutoff, numResOverlap, "0", &verifier);
            
            overlapFinder.insertStructures(batch.first);
            overlapFinder.findOverlaps(batch.second, outFile);
        }

    }
    // Search for overlaps
    /*FuseCandidateFinder fuser(numResOverlap, general, rmsdCutoff, numBatches, batchIndex - 1);
    fuser.minCosAngle = minCosAngle;
    fuser.writeFuseCandidates(binaryPath, outPath);*/

    return 0;
}
