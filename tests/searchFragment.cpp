//
//  quickSearch.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 2/26/19.
//

#include "mstfasst.h"
#include "mstoptions.h"
#include "mstsystem.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Define a fragment from a protein structure, search it against a FASST database, and store the matches");
  op.addOption("p", "The complete protein structure", true);
  op.addOption("sel","The selected residues in the protein structure", true);
  op.addOption("max_matches", "Will search for up to this many matches, given the specified criteria", true);
  op.addOption("config","Path to the configuration file (specifies fasst database and rotamer library",true);
  op.addOption("seq", "If provided, specified positions will be searched with a sequence constraint (should be formatted as a selection)");
  op.addOption("sidechains", "If specified, the sidechains of the matches will be stored along with xthe backbone atoms");
  op.setOptions(argc, argv);
  
  // fasstDB and rotlib are defined at the time of compilation
  configFile config(op.getString("config"));
  
  //// Define query fragment
  //Import structure and selector
  Structure prot(op.getString("p"));
  selector sel(prot); //why is this class name lowercase?
  
  vector<Residue*> query_res = sel.selectRes(op.getString("sel"));
  Structure query(query_res);
  cout << "query residue size " << query.residueSize() << endl;
  
  //make a map that stores the residue* to index-in-structure relationship
  /* note that
   1) this assumes that residue order is preserved between the selection vector and structure
   2) the residue index in the query structure starts with 0
   */
  map<Residue*, int> query_indices;
  for (int i = 0; i < query_res.size(); i++) {
    query_indices.emplace(query_res[i],i);
  }
  
  //// Initialize FASST and search
  cout << "Initializing FASST..." << endl;
  
  FASST F;
  F.setOptions(fasstSearchOptions());
  F.setRedundancyProperty("sim");
  
  cout << "Reading database...";
  auto begin = chrono::high_resolution_clock::now();
  F.readDatabase(config.getDB(),1); //loads backbone atoms and deletes the structure
  auto end = chrono::high_resolution_clock::now();
  cout << "Completed... reading the database took " << chrono::duration_cast<std::chrono::milliseconds>(end-begin).count() << " ms" << endl;
  
  
  // Set query and search parameters
  F.setQuery(query);
  F.setRMSDCutoff(RMSDCalculator::rmsdCutoff(query, 1.1, 15));
  F.setMinNumMatches(0); //there may exist no matches given the specified criteria
  F.setMaxNumMatches(op.getInt("max_matches"));
  
  if (op.isGiven("seq")) {
    cout << "sequence constraint provided..." << endl;
    // add sequence constraints, as needed
    Structure splitQuery = F.getQuery();
    fasstSeqConstSimple seqConst(splitQuery.chainSize());
    
    vector<Residue*> const_residues = sel.selectRes(op.getString("seq"));
    for (int i = 0; i < const_residues.size(); i++) {
      const Residue& res = splitQuery.getResidue(query_indices[const_residues[i]]);
      seqConst.addConstraint(res.getChain()->getIndex(), res.getResidueIndexInChain(), {res.getName()});
    }
    if (seqConst.hasConstraints()) F.options().setSequenceConstraints(seqConst);
  }
  
  // search
  fasstSolutionSet matches;
  matches = F.search();
  
  // Make directory to store the fragments
  string fragDir = "./output/";
  MstSys::cmkdir(fragDir,true);
  
  bool detailed = false;
  if (op.isGiven("sidechains")) {
    cout << "If sidechains exist in the database, they will be included in the output..." << endl;
    detailed = true;
  }
  
  vector<Structure> match_structures;
  F.getMatchStructures(matches,match_structures,detailed);
  
  for (int i = 0; i != match_structures.size(); i++) {
    Structure& s = match_structures[i];
    s.reassignChainsByConnectivity();
    s.writePDB(fragDir + "match" + MstUtils::toString(i) + ".pdb");
  }
  
  return 1;
}
