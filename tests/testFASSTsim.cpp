//
//  testFASSTsim.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 3/24/19.
//

#include "mstsystem.h"
#include "mstoptions.h"
#include "mstfasst.h"

/* import pdb
    initialize FASST objext
 settings
 
 search
 for each match
 rmsd, db_id, position?, sequence? */

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Test the sim property of FASST");
    op.addOption("p", "query structure", true);
    op.setOptions(argc, argv);
    
    Structure query(op.getString("p"));
    
    FASST F;
    fasstSearchOptions options = fasstSearchOptions();
    F.setOptions(options);
    F.setRedundancyProperty("sim");
    
    string fasstDBPath = "/home/grigoryanlab/library/databases/FASST/db-12657.sim";
    cout << "Reading database..." << endl;
    auto begin = chrono::high_resolution_clock::now();
    F.readDatabase(fasstDBPath,1);
    auto end = chrono::high_resolution_clock::now();
    cout << "Reading the database took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
    cout << "Completed..." << endl;
    
    cout << "redundancy cut " << F.getRedundancyCut() << endl;
    cout << "redundancy property " << F.getRedundancyProperty() << endl;
    
    F.setQuery(query);
    F.setRMSDCutoff(1);
    F.setMinNumMatches(20);
    
    for (mstreal sim_cut = 1; sim_cut > 0; sim_cut = sim_cut - .1) {
        options.setRedundancyCut(sim_cut);
        F.setOptions(options);
        F.setQuery(query);
        F.setMinNumMatches(20);
        cout << "redundancy cut " << F.getRedundancyCut() << endl;
        cout << "redundancy property " << F.getRedundancyProperty() << endl;
        
        fasstSolutionSet matches = F.search();
        cout << matches.size() << " matches" << endl;
        
        for (int match_id = 0; match_id != matches.size(); match_id++) {
            fasstSolution match = matches[match_id];
            cout << "rmsd: " << match.getRMSD() << "\t";
            cout << "fasst_id: " << match.getTargetIndex() << "\t";
            cout << endl;
        }
    }
    
    
    
    return 1;
}
