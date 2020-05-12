//
//  exportAtomPositions.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 2/14/20.
//

//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"
#include "msttypes.h"

//tpd dependencies
#include "secondarystructure.h"
#include "utilities.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Load seeds from a binary file and write their coordinates, plus other information to a CSV");
  op.addOption("extfrag_bin", "If a binary file of extended fragment structures already exists, will skip creating fragments/extending them and just compare the provided ones to the peptide.",true);
  op.addOption("out", "The prefix for the csv file which will contain the CA coordinates/coloring",true);
  op.addOption("sample", "The fraction of seeds that should be exported (selected randomly)");
  op.setOptions(argc, argv);
  
  
  // Variables provided by user
  string extfrag_bin = op.getString("extfrag_bin");
  string out_prefix = op.getString("out");
  
  fstream out;
  string filePath = out_prefix + "_seed_ca.lcloud";
  MstUtils::openFile(out, filePath, fstream::out);
  
  secondaryStructureClassifier classifier;
  classifier.writeCaInfotoLineFile(extfrag_bin, 100000, out);
  
  out.close();
  
  return 0;
}
