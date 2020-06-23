//
//  testBinMetaData.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/22/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include "structure_iter.h"
#include "mstoptions.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Read the structures/meta data in a version 1 seed binary file");
  op.addOption("bin", "path to binary file",true);
  op.addOption("top", "only look at top N structures",true);
  op.setOptions(argc,argv);
  
  StructuresBinaryFile seeds(op.getString("bin"),true,1);
  seeds.scanFilePositions();
  seeds.reset();
  
  int top = op.getInt("top",10000000);
  int count = 0;
  while (seeds.hasNext() && count < top) {

    Structure* extended_fragment = seeds.next();
    
    cout << "name: " << extended_fragment->getName() << endl;
    cout << "seq: " << seeds.getStructurePropertyInt("seq",extended_fragment->getName()) << endl;
    cout << "match_rmsd: " << seeds.getStructurePropertyReal("match_rmsd",extended_fragment->getName()) << endl;
    cout << "rmsd_adjust: " << seeds.getStructurePropertyReal("rmsd_adj",extended_fragment->getName()) << endl;
    
    delete extended_fragment;
    
    count += 1;
  }
  return 1;
}
