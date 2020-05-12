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


int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Loads a binary file, finds the structure with the provided name, and saves it as a PDB file");
  op.addOption("name", "Name of the structure in the binary file");
  op.addOption("list", "The path to a file listing names of structures in the binary file");
  op.addOption("bin", "The path to the binary file with the structure", true);
  op.setOptions(argc, argv);
  
  if ((!op.isGiven("name")) && (!op.isGiven("list"))) MstUtils::error("--name or --list must be provided");
  
  string bin_path = op.getString("bin");
  
  vector<string> structure_names;
  if (op.isGiven("name")) {
    structure_names.push_back(op.getString("name"));
  } else {
    structure_names = MstUtils::fileToArray(op.getString("list"));
  }
  
  StructuresBinaryFile bin(bin_path);
  
  for (string name : structure_names) {
    Structure* S = bin.getStructureNamed(name);
    S->writePDB(name+".pdb");
    if (S != NULL) delete S;
  }
}
