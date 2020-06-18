//
//  seedCentroidDistance.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/13/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//


#include "mstsystem.h"
#include "mstoptions.h"

#include "structure_iter.h"
#include "utilities.h"
#include "benchmarkutilities.h"

//using namespace std;
//using namespace std::chrono;

int main (int argc, char *argv[]) {
  
  // Get command-line arguments
  MstOptions op;
  op.setTitle("Takes a set of complexes and seeds generated around them and finds the average distance distribution");
  op.addOption("list", "path to a file where each line is '/path/to/seed_bin /path/to/complex'", true);
  op.setOptions(argc, argv);
  
  cout << "test" << endl;
  
  mstreal min_value = 0.0;
  mstreal max_value = 50.0;
  int num_bins = 150.0;
  
  string hist_file = "seed_distance_histogram.csv";
  
  string list = op.getString("list");
  
  seedCentroidDistance seedDist(list, min_value, max_value, num_bins);
  
  histogram hist = seedDist.getHist();
  
  hist.writeHistFile(hist_file);
  
  return 0;
}
