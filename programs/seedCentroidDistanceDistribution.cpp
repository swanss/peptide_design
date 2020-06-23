//
//  seedCentroidDistanceDistribution.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/19/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include "mstoptions.h"
#include "mstoptim.h"
#include "msttypes.h"
#include "benchmarkutilities.h"

int main(int argc, char* argv[]) {
  MstOptions op;
  op.setTitle("Samples seeds from multiple binary files, computes their distance from the protein, and averages these to build a histogram");
  op.addOption("list","a file where each line is the path to a seed binary file");
  op.setOptions(argc, argv);
  
  string list = op.getString("list");
  mstreal min_val = 0.0;
  mstreal max_val = 50.0;
  int num_bins = 150;
  
  seedCentroidDistance dist(list,min_val,max_val,num_bins);
  
  histogram hist = dist.getHist();
  
  hist.writeHistFile("seed_centroid_distance.csv");

  return 1;
}
