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
    opts.addOption("out", "Directory into which to write overlaps", true);
    opts.addOption("overlapSize", "Number of residues that must overlap between two residues (default 2)", false);
    opts.addOption("overlapRMSD", "Maximum RMSD between sets of residues to be considered overlapping (default 1.0)", false);
    opts.addOption("minCosAngle", "Minimum cosine angle between Ca vectors allowed in a potential overlap (default -1.0)", false);
    opts.addOption("batch", "The batch number (from 1 to numBatches)", false);
    opts.addOption("numBatches", "The number of batches", false);
    opts.setOptions(argc, argv);
    
    //string filesPath = opts.getString("files");
    //string dataPath = opts.getString("data");
    string binaryPath = opts.getString("bin");
    string outPath = opts.getString("out");
    if (!MstSys::fileExists(outPath))
        MstSys::cmkdir(outPath);

    int numResOverlap = opts.getInt("overlapSize", 2);
    float rmsdCutoff = opts.getReal("overlapRMSD", 1.0);
    float minCosAngle = opts.getReal("minCosAngle", -1.0);

    int batchIndex = opts.getInt("batch", 1);
    int numBatches = opts.getInt("numBatches", 1);
    if (batchIndex < 1 || batchIndex > numBatches) {
        cerr << "Batch index must be between 1 and numBatches" << endl;
        return 1;
    }

    // Search for overlaps
    /*SeedListFile seedFile(filesPath);
    auto fileContents = seedFile.read(dataPath);
    vector<string> candidatePaths = fileContents.first;*/
    FuseCandidateFinder fuser(numResOverlap, general, rmsdCutoff, numBatches, batchIndex - 1);
    fuser.minCosAngle = minCosAngle;
    // Seed chain is always "0"
    //vector<string> chainIDs(candidatePaths.size(), "0");
    //fuser.writeFuseCandidates(candidatePaths, &chainIDs, outPath, dataPath);
    fuser.writeFuseCandidates(binaryPath, outPath);

    return 0;
}
