//
//  testFASSTDBProperties.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 10/10/19.
//

#include "msttypes.h"
#include "mstfasst.h"
#include "mstoptions.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Reports which properties are defined in a FASSTDB");
  op.addOption("dbPath","The path to the DB",true);
  op.addOption("mem_save","0: no memory save, 1: no sidechains, 2: destroy structures after getting atoms",true);
  op.setOptions(argc, argv);
  
  int mem_save = op.getInt("mem_save");
  
  //import DB
  FASST F;
  F.readDatabase(op.getString("dbPath"),mem_save);
  
  //Structure properties
  int target_num = F.numTargets();
  int chain_num = 0;
  for (int target_id = 0; target_id < target_num; target_id++) {
    Structure* S = F.getTarget(target_id);
    chain_num += S->chainSize();
  }
  
  cout << "Structures in DB: " << target_num << endl;
  cout << "Chains in DB: " << chain_num << endl;
  
  cout << "Residue property: phi = " << F.isResiduePropertyDefined("phi") << endl;
  cout << "Residue property: psi = " << F.isResiduePropertyDefined("psi") << endl;
  cout << "Residue property: omega = " << F.isResiduePropertyDefined("omega") << endl;
  cout << "Residue property: env = " << F.isResiduePropertyDefined("env") << endl;
  cout << "Residue pair property: cont = " << F.isResiduePairPropertyPopulated("cont") << endl;
  cout << "Residue pair property: interfering = " << F.isResiduePairPropertyPopulated("interfering") << endl;
  cout << "Residue pair property: interfered = " << F.isResiduePairPropertyPopulated("interfered") << endl;
  cout << "Residue pair property: bb = " << F.isResiduePairPropertyPopulated("bb") << endl;
  cout << "Residue relationship property sim: " << MstUtils::toString(F.isResidueRelationshipPopulated("sim")) << endl;
  
  return 1;
}
