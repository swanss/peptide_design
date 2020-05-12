//
//  testSecStructure.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 2/17/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"
#include "msttypes.h"

//tpd dependencies
#include "secondarystructure.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Load PDB, classify secondary structure of each residues, and write to file");
  op.addOption("pdb", "Structure of interest",true);
  op.addOption("out", "The prefix for the sequence file. The suffix is always '_sec_struct_out.tsv'",true);
  op.setOptions(argc, argv);
  
  // Variables provided by user
  string pdb = op.getString("pdb");
  string out_prefix = op.getString("out");
  
  Structure S(pdb);
  
  fstream out;
  string filePath = out_prefix + "_sec_struct_out.tsv";
  MstUtils::openFile(out, filePath, fstream::out);
  
  for (int flank = 0; flank <= 10; flank++) {
    int max_agreement = (flank*2)+1; //max is all flanking residues required to agree
    for (int agreement = 0; agreement < max_agreement; agreement++) {
      secondaryStructureClassifier classifier(flank, agreement);
      classifier.set_stride(true);
      out << flank << "\t" << agreement << "\t" << classifier.classifyStruct(S) << endl;
    }
  }
  out.close();
  
  return 0;
}
