//
//  testOverlapPerformance.cpp
//
//  Created by Venkatesh Sivaraman on 7/1/20.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "overlaps.h"
#include "seedgraph.h"
#include "seedscore.h"

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Tests the performance of overlap finding at different seed densities.");
    opts.addOption("bin", "Path to a binary file containing seed structures", true);
    opts.addOption("id", "Identifier for this test run (used to print results)", true);
    opts.addOption("overlapSize", "Number of residues that must overlap between two residues (default 2)", false);
    opts.addOption("deviation", "Maximum deviation allowed between individual residues in an overlap segment (default 1.0)", false);
    opts.addOption("minCosAngle", "Minimum cosine angle between residue normal vectors in an overlap segment (default -1.0)", false);
    opts.addOption("cap", "Maximum expected number of seeds to use in a single run (default 400,000)", false);
    opts.setOptions(argc, argv);
    
    string binaryPath = opts.getString("bin");
    string outPath = opts.getString("out");
    StructuresBinaryFile binaryFile(binaryPath);

    int numResOverlap = opts.getInt("overlapSize", 2);
    float maxDeviation = opts.getReal("deviation", 1.0);
    float minCosAngle = opts.getReal("minCosAngle", -1.0);
    
    OverlapVerifier *verifier = new MaxDeviationVerifier(maxDeviation);
    if (minCosAngle > -1.0) {
        verifier = new CompositeVerifier(verifier, new NormalVectorVerifier(minCosAngle));
    }

    vector<int> moduli = { 20, 15, 10, 5, 4, 3, 2, 1 };
    
    vector<int> numOverlaps;
    vector<double> hashTimes;
    vector<double> searchTimes;
    vector<int> numSeeds;
    vector<int> numResidues;
    vector<double> volumes;
    
    int totalStructureCount = binaryFile.structureCount();
    int structureCap = opts.getInt("cap", 400000);
    
    for (int modulus: moduli) {
        if (totalStructureCount / modulus >= structureCap) {
            cout << "Modulus " << modulus << " would take too long, stopping" << endl;
            break;
        }
        
        cout << "Modulus: " << modulus << endl;
        
        mstreal xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9, zlo = 1e9, zhi = -1e9;
        
        binaryFile.reset();
        vector<Structure *> structures;
        int residueCount = 0;
        while (binaryFile.hasNext()) {
            if (rand() % modulus != 0) {
                binaryFile.skip();
                continue;
            }
            
            Structure *seed = binaryFile.next();
            structures.push_back(seed);
            
            if (seed->residueSize() == 0)
                continue;
            
            Structure floatingChains;
            extractChains(*seed, "0", floatingChains);
            
            mstreal ixlo, ixhi, iylo, iyhi, izlo, izhi;
            ProximitySearch::calculateExtent(floatingChains, ixlo, iylo, izlo, ixhi, iyhi, izhi);
            xlo = min(xlo, ixlo);
            xhi = max(xhi, ixhi);
            ylo = min(ylo, iylo);
            yhi = max(yhi, iyhi);
            zlo = min(zlo, izlo);
            zhi = max(zhi, izhi);
            
            residueCount += floatingChains.residueSize();
        }
        
        numSeeds.push_back(structures.size());
        numResidues.push_back(residueCount);
        
        vector<mstreal> bbox = { xlo, xhi, ylo, yhi, zlo, zhi };
        volumes.push_back((xhi - xlo) * (yhi - ylo) * (zhi - zlo));

        // Measure time to compute overlaps
        MstTimer timer;
        timer.start();

        CAResidueHasher<> hasher(bbox, 0.3);
        OverlapFinder<CAResidueHasher<>> overlapFinder(hasher, maxDeviation, numResOverlap, "0", verifier);
        overlapFinder.verbose = false;
        
        overlapFinder.insertStructures(structures);
        timer.stop();
        double performance = ((double)timer.getDuration(MstTimer::msec) / 1000.0);
        hashTimes.push_back(performance);
        cout << "Hash time: " << performance << endl;
        
        timer.start();
        vector<FuseCandidate> results = overlapFinder.findOverlaps(structures, true);
        timer.stop();
        cout << "Num overlaps: " << results.size() << endl;
        numOverlaps.push_back(results.size());
        
        performance = ((double)timer.getDuration(MstTimer::msec) / 1000.0);
        
        cout << "Search time: " << performance << endl;
        searchTimes.push_back(performance);
        
        for (Structure *s: structures)
            delete s;
    }
    
    cout << "-\tID\tModulus\tNum Seeds\tNumResidues\tBBox Volume\tNum Overlaps\tHash Time\tSearch Time" << endl;
    for (int i = 0; i < searchTimes.size(); i++) {
        cout << "RESULTS\t" << opts.getString("id") << "\t" << moduli[i] << "\t" << numSeeds[i] << "\t" << numResidues[i] << "\t" << volumes[i] << "\t" << numOverlaps[i] << "\t" << hashTimes[i] << "\t" << searchTimes[i] << endl;
    }
    
    delete verifier;
    
    return 0;
}
