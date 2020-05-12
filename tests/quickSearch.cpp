//
//  searchFragment.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 2/26/19.
//

// searches query and stores 10 matches


#include "mstsystem.h"
#include "mstoptions.h"
#include "mstfasst.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Search query against database and store matches");
  op.addOption("p", "query structure", true);
  op.addOption("config","Path to the configuration file (specifies fasst database and rotamer library",true);
  op.setOptions(argc, argv);
  
  string fragDir = "./output/fragments/";
  MstSys::cmkdir(fragDir,true);
  
  // Import structure
  Structure target(op.getString("p"));
  
  configFile config(op.getString("config"));
  
  // Initialize fragmenter object
  cout << "Initializing FASST object..." << endl;
  
  FASST F;
  fasstSearchOptions foptsBase;
  fasstSolutionSet matches;
  
  F.setOptions(fasstSearchOptions());
  foptsBase = F.options();
  //    F.setRedundancyProperty("sim");
  cout << "Reading database..." << endl;
  auto begin = chrono::high_resolution_clock::now();
  F.readDatabase(config.getDB(),2);
  auto end = chrono::high_resolution_clock::now();
  cout << "Reading the database took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  cout << "Completed..." << endl;
  
  F.setQuery(target);
  F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(target, 1.1, 15));
  F.setMinNumMatches(1);
  F.setMaxNumMatches(10);
  matches = F.search();
  
  for (int i = 0; i != matches.size(); i++) {
    fasstSolution sol;
    sol = matches[i];
    Structure match_protein;
    match_protein = F.getMatchStructure(sol,false,FASST::matchType::WITHGAPS,true);
    match_protein.writePDB("match" + MstUtils::toString(i) + ".pdb");
    cout << "print residues: ";
    for (int i = 0; i != match_protein.residueSize(); i++) {
      cout << match_protein.getResidue(i) << " ";
    }
    cout << endl;
    
  }
  
  return 1;
}
