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
    mstreal max_val = 25.0;
    int num_bins = 100;
    int sample_n = 1000000;
    
    //get the list of binary files/structures
    //e.g. path/to/binary path/to/structure
    vector<string> files = MstUtils::fileToArray(list);
    cout << "there are " << files.size() << " sets of seeds" << endl;
    
    if (files.size() * sample_n > INT_MAX) MstUtils::error("Note: it is possible that the counts in a single bin could exceed the maximum possible value of an 'int'. If this occurs, consider sampling less seeds or refactoring");
    
    vector<histogram> all_histograms;
    for (string file_name : files) {
        cout << file_name << endl;
        vector<string> split = MstUtils::split(file_name);
        string bin_path = split[0];
        string struct_path = split[1];
        
        //get the peptide chain name
        split = MstUtils::split(struct_path,"_");
        string p_id = split[split.size() - 3]; //should get the peptide chain id
        
        cout << "pdb path: " << struct_path << " bin path: " << bin_path << " peptide chain ID: " << p_id << endl;
        
        Structure complex(struct_path);
        seedStatistics seedStat(complex,p_id,bin_path);
        
        all_histograms.push_back(seedStat.generateDistanceHistogram(min_val, max_val, num_bins, sample_n));
    }
    
    //average the histograms
    histogram summary_hist(min_val,max_val,num_bins);
    for (int i = 0; i < num_bins; i++) {
        mstreal bin_sum = 0;
        for (int j = 0; j < all_histograms.size(); j++) bin_sum += all_histograms[j].getVal(i);
        mstreal bin_mean = bin_sum / all_histograms.size();
        summary_hist.setBinVal(i,bin_mean);
    }
    
    summary_hist.writeHistFile("seed_centroid_distance_normalized.csv");
    
    return 1;
}
