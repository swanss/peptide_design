//
//  testFASSTPerformance.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 3/4/19.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "fusehelper.h"
#include "seedgraph.h"
#include "seedscore.h"

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    string dataPath = "/home/ifsdata/scratch/grigoryanlab/swans/peptide_design_benchmark/seed_retrieval/target_centered_approach/retrieveSeedsTest3/1DUY/output/seeds/";
    string fasstDB = "/home/grigoryanlab/library/databases/dTERMen.databases/2019-01-22/dtermen.sim";
    Structure *target = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb");
    Structure *seed = new Structure(dataPath + "1DUY_D26/1DUY_D26_match_158_seed_0.pdb");
    
    string *queryWritePath = nullptr;
    if (argc > 1) {
        queryWritePath = new string(argv[1]);
        cout << "query write path: " << *queryWritePath << endl;
    }

    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    
    SequenceCompatibilityScorer scorer(target, rParams, cParams, fasstDB, 2, 0, 0.4, 0.05, 0.25, 1, 8000, 0.7);
    scorer.queryWritePath = queryWritePath;
    cout << "Loaded scorer" << endl;
    high_resolution_clock::time_point startTime = high_resolution_clock::now();
    auto result = scorer.score(seed);
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    cout << "Scoring took " << duration_cast<seconds>( endTime - startTime ).count() << " seconds" << endl;
    cout << "Result: ";
    for (auto it: result) {
        cout << it.first->getName() << it.first->getResidueIndex() << " = " << it.second << ", " << endl;
    }
    cout << endl;
    
    if (queryWritePath != nullptr)
        delete queryWritePath;

    return 0;
}
