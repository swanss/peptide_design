//
//  getFromBin.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 4/4/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//
//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "utilities.h" //for generateAllCombinationsKRes
#include "structure_iter.h"


int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Loads a binary file, adds the structure to the file.");
  op.addOption("pdb", "path to the structure to add to the binary file",true);
  op.addOption("name", "name of the structure",true);
  op.addOption("bin", "The path to the binary file with the structure", true);
  op.setOptions(argc, argv);
  
  string bin_path = op.getString("bin");
  string pdb = op.getString("pdb");
  string name = op.getString("name");
  
  Structure* S = new Structure(pdb);
  S->setName(name);
  
  //structure can only have one chain (it will be interpreted as a seed chain)
  MstUtils::assert(S->chainSize() == 1,"There should only by a single chain in the input file");
  Chain* C = &S->getChain(0);
  C->setID("0");
  
  //open in write mode to add the structure
  bool read = false;
  bool append = true;
  StructuresBinaryFile* bin = new StructuresBinaryFile(bin_path,read,append);
  
  bin->appendStructure(S);
  
  delete S;
  delete bin;
  
  //open in read mode to check that the structure was added
  bin = new StructuresBinaryFile(bin_path,true);
  
  S = bin->getStructureNamed(name);
  cout << "Structure with name: " << S->getName() << " now in the binary file" << endl;
  
  return 0;
}
