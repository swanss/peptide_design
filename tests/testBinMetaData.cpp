//
//  testBinMetaData.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/22/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include "structure_iter.h"

int main() {
  StructuresBinaryFile seeds("/Users/sebastianswanson/Keating/tests/1LB5_benchmark/termextension_output/extendedfragments.bin",true,1);
  seeds.scanFilePositions();
  seeds.reset();
  
  
  while (seeds.hasNext()) {

    Structure* extended_fragment = seeds.next();
    
    cout << "name: " << extended_fragment->getName() << endl;
    cout << "seq: " << seeds.getStructurePropertyInt("seq",extended_fragment->getName()) << endl;
    cout << "match_rmsd: " << seeds.getStructurePropertyReal("match_rmsd",extended_fragment->getName()) << endl;
    cout << "rmsd_adjust: " << seeds.getStructurePropertyReal("rmsd_adj",extended_fragment->getName()) << endl;
    
    delete extended_fragment;
  }
  return 1;
}
