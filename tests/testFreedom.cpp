//
//  testFreedom.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 3/1/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include "mstsystem.h"
#include "mstoptions.h"
#include "mstcondeg.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("For visualizing the results of the freedom calculation in ConFind.");
  op.addOption("p", "The structure which will be used to compute contact information",true);
  op.setOptions(argc,argv);
  
  string RL = "/Users/sebastianswanson/Keating/utilities/repos/MST/testfiles/rotlib.bin";
  //  string RL = "/home/ifs-users/swans/MST_workspace/MST/testfiles/rotlib.bin";
  
  Structure S(op.getString("p"));
  ConFind C(RL,S,false);
  
  vector<Residue*> all_residues = S.getResidues();
  
  for (Residue* R : all_residues) {
    mstreal freedom = C.getFreedom(R);
    cout << *R << " freedom: " << freedom << endl;
    for (Atom* A : R->getAtoms()) A->setB(freedom);
  }
  
  S.writePDB("test_freedom.pdb");
}

