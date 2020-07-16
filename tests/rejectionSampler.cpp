//
//  rejectionSampler.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/29/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include "mstoptions.h"
#include "msttypes.h"

#include "benchmarkutilities.h"

int main(int argc, char* argv[]) {
  MstOptions op;
  op.setTitle("Test program to verify that rejection sampling is working properly");
  op.addOption("hist", "Path to a file that can be read by Histogram");
  op.setOptions(argc, argv);
  
  string hist_path = op.getString("hist");
  
  rejectionSampler rSampler(hist_path);
  
  fstream out;
  MstUtils::openFile(out,"out.tsv",ios_base::out);
  out << "number\tvalue" << endl;
  int numFinalSampled = 1000; int count = 0; mstreal val;
  for (int i = 0; i < numFinalSampled; i++)  {
    bool sampled = false;
    while (!sampled) {
      //sample value
      val = MstUtils::randUnit() * 25;
      cout << "sample: " << val << endl;
      //check if rejected
      sampled = rSampler.accept(val);
    }
    cout << "accept: " << val << endl;
    out << i << "\t";
    out << val << endl;
  }

  cout << "Done" << endl;
  return 1;
}
