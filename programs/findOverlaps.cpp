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

using namespace std;

int main (int argc, char *argv[]) {
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Finds overlaps between all pairs of seeds.");
    opts.addOption("seedBin", "Path to a binary file containing seed structures.", true);
    opts.addOption("out", "Path to file at which to write overlaps ('.csv prefix will be added automatically')'", true);
    opts.addOption("overlapSize", "Number of residues that must overlap between two residues. Must be even (default 4)", false);
    opts.addOption("maxDeviation", "Distance cutoff (Å) between alpha-carbons in an overlap segment (default 1.0)", false);
    opts.addOption("minCosAngle", "Cosine angle threshold between residue normal vectors in an overlap segment (default 0.5)", false);
    opts.addOption("worker", "The index of this worker, from 1 to numWorkers. (default: 1)", false);
    opts.addOption("numWorkers", "The number of workers (default: 1)", false);
    opts.addOption("batchSize", "The number of structures to use in each batch. (default 200,000)", false);
    opts.addOption("bruteForce", "If provided, use a brute-force all-to-all comparison", false);
    opts.addOption("limitBatches", "Number of batches to run (to produce a smaller debugging set)", false);
    opts.addOption("mock", "If provided, do not perform any overlap calculations, just list the batches that would be computed", false);
    opts.setOptions(argc, argv);
    
    //string filesPath = opts.getString("files");
    //string dataPath = opts.getString("data");
    string binaryPath = opts.getString("seedBin");
    string outPath = opts.getString("out") + ".csv";

    int numResOverlap = opts.getInt("overlapSize", 4);
    float maxDeviation = opts.getReal("maxDeviation", 1.0);
    float minCosAngle = opts.getReal("minCosAngle", 0.5);

    int worker = opts.getInt("worker", 1);
    int numWorkers = opts.getInt("numWorkers", 1);
    int batchSize = opts.getInt("batchSize", 200000);
    int limitBatches = opts.getInt("limitBatches", -1);
    
    if (numResOverlap % 2 != 0) cout << "Warning: overlapSize is not even and can't be used to generate a seed graph" << endl;
    
    if (worker < 1 || worker > numWorkers) {
        cerr << "Batch index must be between 1 and numWorkers" << endl;
        return 1;
    }
    cout << "Worker " << worker << " of " << numWorkers << " total workers"<< endl;
    
    OverlapVerifier *verifier = new MaxDeviationVerifier(maxDeviation);
    if (minCosAngle > -1.0) {
        verifier = new CompositeVerifier(verifier, new NormalVectorVerifier(minCosAngle));
    }

    FuseCandidateFile outFile(outPath);
    BatchPairStructureIterator structureIter(binaryPath, worker - 1, numWorkers, batchSize); // large batch size
    
    // Measure time to compute overlaps
    MstTimer timer;
    timer.start();
    
    int batchesRun = 0;
    
    while (structureIter.hasNext()) {
        if (limitBatches >= 0 && batchesRun > limitBatches)
            break;
        
        auto batch = structureIter.next();
        cout << batch.first.size() << " " << batch.second.size() << endl;
        
        if (opts.isGiven("mock")) {
            cout << "Batch: " << batch.first[0]->getName() << ", " << batch.second[0]->getName() << endl;
            batchesRun++;
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
                        if (s2 == s1) continue;
                        
                        Chain *c2 = s2->getChainByID("0");
                        if (!c2)
                            continue;
                        
                        vector<Residue *> residues2 = c2->getResidues();
                        for (int resIdx2 = 0; resIdx2 < residues2.size() - numResOverlap + 1; resIdx2++) {
                            vector<Residue *> segment2(residues2.begin() + resIdx2, residues2.begin() + resIdx2 + numResOverlap);
                            
                            if (verifier->verify(segment1, segment2)) {
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
                batchesRun++;
            }
        } else {
            mstreal xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9, zlo = 1e9, zhi = -1e9;
            vector<Structure *> &firstStructures = batch.first;
            vector<Structure *> &secondStructures = batch.second;
            
            auto it = firstStructures.begin();
            while (it != secondStructures.end()) {
                Structure *seed = *it;
                if (seed->residueSize() == 0) {
                      cout << "seed has no residues" << endl;
                      continue;
                }
                
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
            OverlapFinder<CAResidueHasher<>> overlapFinder(hasher, maxDeviation, numResOverlap, "0", verifier);
            
            overlapFinder.insertStructures(batch.first);
            // Tell the overlap finder whether to avoid double-counting symmetric overlaps
            bool symmetricBatch = firstStructures[0]->getName() == secondStructures[0]->getName();
            if (symmetricBatch) cout << "Symmetric batch" << endl;

            overlapFinder.findOverlaps(batch.second, outFile, symmetricBatch);
            batchesRun++;
        }
    }
    if (batchesRun == 0) MstUtils::error("No batches assigned to worker with index: "+MstUtils::toString(worker),"findOverlaps::main");
    
    // Search for overlaps
    /*FuseCandidateFinder fuser(numResOverlap, general, maxDeviation, numWorkers, worker - 1);
    fuser.minCosAngle = minCosAngle;
    fuser.writeFuseCandidates(binaryPath, outPath);*/

    timer.stop();
    cout << "Elapsed time: " << ((double)timer.getDuration(MstTimer::msec) / 1000.0) << endl;

    delete verifier;
    
    cout << "Done!" << endl;
    return 0;
}
