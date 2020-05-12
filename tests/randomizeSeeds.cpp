//
//  randomizeSeeds.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 3/5/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "coverage.h"
#include "secondarystructure.h"
#include "seedutilities.h"
#include "utilities.h"

////structgen dependencies
#include "Util.h"
#include "vdwRadii.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Load seeds from a binary file, randomize their position and orientation, and write them to a new binary file. Next, determine the coverage.");
  op.addOption("pdb", "The peptide-protein complex that the seeds cover",true);
  op.addOption("p_id", "The ID of the peptide chain",true);
  op.addOption("bin_in", "Binary file with the current seeds",true);
  op.addOption("out", "The prefix for the new binary/csv files",true);
  op.addOption("config", "The path to the config file",true);
  op.addOption("position", "If provided, will randomize the position of seeds");
  op.addOption("orientation","If provided, will randomize the orientation of seeds");
  op.addOption("sample","If provided, will sample an equal length segment from the PDB for each seed");
  op.setOptions(argc, argv);
  
  int seed_num = 100000;
  
  // Variables provided by user
  string bin_in = op.getString("bin_in");
  string out_prefix = op.getString("out");
  string pdb_path = op.getString("pdb");
  string p_id = op.getString("p_id");
  bool sample = op.isGiven("sample");
  bool position = op.isGiven("position");
  bool orientation = op.isGiven("orientation");
  
  configFile config(op.getString("config"));
  string RL_path = config.getRL();
  string DB_path = config.getDB();
  
  //// Make directories
  bool makeParents = true;
  // Fragment output folder
  string outDir = "termextension_output/";
  string covDir = "coverage_output/";
  MstSys::cmkdir(outDir,makeParents);
  MstSys::cmkdir(covDir,makeParents);
  
  //names
  string bin_out = outDir + out_prefix + ".bin";
  string lcloud_out = outDir + out_prefix + "_random_seed_ca.lcloud";
  
  
  Structure S(pdb_path);
  
  MstTimer timer;
  timer.start();
  if (sample) {
    cout << "seeds from binary file" << endl;
    naiveSeedsFromBin naiveSeeds(S,p_id);
    naiveSeeds.newPose(bin_in,outDir,out_prefix,position,orientation);
  } else {
    cout << "seeds from DB" << endl;
    naiveSeedsFromDB naiveSeeds(S,p_id,DB_path);
    naiveSeeds.newPose(bin_in,outDir,out_prefix,position,orientation);
  }
  timer.stop();
  
  cout << "Randomized all seeds in " << timer.getDuration() << " seconds" << endl;
  
  // write the seed CA atoms to a file for visualization
  fstream out;
  MstUtils::openFile(out, lcloud_out, fstream::out);
    
  secondaryStructureClassifier classifier;
  classifier.writeCaInfotoLineFile(bin_out, seed_num, out);

  out.close();
  
  //get the coverage
  mstreal max_rmsd = 3.0;
  int max_seed_length = 10;
  interfaceCoverage IC(S, p_id, max_rmsd, max_seed_length, RL_path);
  IC.findCoveringSeeds(bin_out);
  IC.writeCoverageToFiles(covDir);
  
  //report statistics
  seedStatistics stats(S, p_id);
  stats.writeStatisticstoFile(bin_out, outDir, out_prefix, seed_num);
  
  return 0;
}

