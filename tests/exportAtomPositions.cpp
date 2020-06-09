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
  op.addOption("list", "Path to a file where each line specifies a seed in the binary file that should be exported");
  op.addOption("sample", "The number of seeds that should be sampled from the binary file");
  op.setOptions(argc, argv);
  
  if (!op.isGiven("list") || !op.isGiven("sample")) MstUtils::error("Either '--list' or '--sample' must be provided (but not both)");
  if (op.isGiven("list") && op.isGiven("sample")) MstUtils::error("'--list' and '--sample' cannot both be provided");
  
  // Variables provided by user
  string extfrag_bin = op.getString("extfrag_bin");
  string out_prefix = op.getString("out");
  
  fstream out;
  string filePath = out_prefix + "_seed_ca.lcloud";
  MstUtils::openFile(out, filePath, fstream::out);
  
  secondaryStructureClassifier classifier;
  if (op.isGiven("sample")){
    int sample = op.getInt("sample",100000);
    classifier.writeCaInfotoLineFile(extfrag_bin, sample, out);
  } else {
    string list_file = op.getString("list");
    vector<string> seed_names = MstUtils::fileToArray(list_file);
    classifier.writeCaInfotoLineFile(extfrag_bin, seed_names, out);
  };
  
  out.close();
  return 0;
};
