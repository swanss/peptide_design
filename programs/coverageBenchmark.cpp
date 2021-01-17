//
//  coverageBenchmark.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/27/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "utilities.h"
#include "coverage.h"
#include "benchmarkutilities.h"
#include "termextension.h"
#include "secondarystructure.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Maps seeds to the peptide in a peptide-protein complex and computes various statistics.");
    op.addOption("bin_path", "path to a seed binary file.",true);
    op.addOption("pdb", "Structure file containing peptide chain and protein chain(s)", true);
    op.addOption("peptide", "Peptide chain ID", true);
    op.addOption("max_rmsd", "The max RMSD threshold used when determining whether a seed aligns to the peptide or not. (default 2.0)");
    op.addOption("config","Path to the configuration file (specifies fasst database and rotamer library)",true);
    op.addOption("hist","Path to the histogram file with the seed distance distribution that will be matched in the null model seeds");
    op.addOption("no_clash_check", "Filter seed poses by clashes to the protein");
    op.setOptions(argc, argv);
    
    MstTimer timer;
    
    // Variables set at the time of compilation
    int max_seed_length = 10; //this limits the length of kmers that are compared between the peptide and a seed
    int num_sampled = 10000; //this is approximately the number of seeds that will be sampled when writing line clouds
    
    // Variables provided by user
    string extfrag_bin = op.getString("bin_path");
    Structure complex(op.getString("pdb"));
    string p_cid = op.getString("peptide");
    mstreal max_rmsd = op.getReal("max_rmsd",2.0);
    configFile config(op.getString("config"));
    string hist_path = op.getString("hist","");
    
    // Set the sequence of the peptide to "unknown"
    selector sel(complex);
    vector<Residue*> peptide_res = sel.selectRes("chain "+p_cid);
    cout << "Selected " << peptide_res.size() << " peptide residues to be renamed to 'UNK'" << endl;
    for (Residue* R : peptide_res) {
        R->setName("UNK");
    }
    
    // Make directories
    bool makeParents = true;
    // Fragment output folder
    string outDir = "output/";
    string covDir = "coverage_output/";
    MstSys::cmkdir(outDir,makeParents);
    MstSys::cmkdir(covDir,makeParents);
    
    //Initialize the interface coverage class
    interfaceCoverage IC(complex, p_cid, config.getRL());
    IC.setMaxRMSD(max_rmsd);
    IC.setMaxSegmentLength(max_seed_length);
    
    //Randomize seeds
    /*
     Two types of randomized seeds
     Type 1 randomized orientation, position sampled from previous seed
     Type 2 randomized position and orientation
     
     To generate type 2, there are two steps:
     1) Generate proposal distribution: sample seeds from the DB and record centroid distance distribution after clash filtering
     2) Sample seeds using rejection sampling: repeat step 1, but accept/reject s.t. the target distribution is recapitulated.
     */
    
    string type1_name = "type1_seeds";
    string type1_bin = outDir + type1_name + ".bin";
    string type2_proposal_name = "type2_proposal_seeds"; //note: the proposal distribution is generated without rejection sampling
    string type2_proposal_bin = outDir + type2_proposal_name + ".bin";
    string type2_name = "type2_seeds";
    string type2_bin = outDir + type2_name + ".bin";
    seedStatistics stats(complex, p_cid);
    if (hist_path != "") {
        bool position = false;
        bool orientation = true;
        
        cout << "Generating null model seeds..." << endl;
        
        //generate type 1 seeds
        naiveSeedsFromDB naiveSeeds(complex, p_cid, extfrag_bin, config.getDB(), config.getRL());
        if (op.isGiven("no_clash_check")) naiveSeeds.setClashChecking(false);
        
        timer.start();
        naiveSeeds.newPose(outDir, type1_name, position, orientation);
        timer.stop();
        cout << "Generated type 1 seeds in " << timer.getDuration() << " seconds" << endl;
        
        //generate type 2 seeds and find distance distribution
        position = true;
        int num_seeds = 1000000; //need to sample this many seeds to get a smooth distribution
        timer.start();
        naiveSeeds.newPose(outDir, type2_proposal_name, position, orientation, num_seeds);
        timer.stop();
        cout << "Generated type 2 seeds (to find proposal distribution) in " << timer.getDuration() << " seconds" << endl;
        
        //generate a proposal histogram
        stats.setBinaryFile(type2_proposal_bin);
        stats.writeStatisticstoFile(outDir, type2_proposal_name, num_sampled);
        histogram proposal = stats.generateDistanceHistogram();
        string proposal_hist_path = outDir + "seed_centroid_distance_proposal.csv";
        proposal.writeHistFile(proposal_hist_path);
        
        //generate type 2 null model seeds with rejection sampling
        histogram target(hist_path);
        rejectionSampler rSampler(proposal,target);
        naiveSeeds.setRejectionSampler(&rSampler);
        timer.start();
        naiveSeeds.newPose(outDir, type2_name, position, orientation);
        timer.stop();
        cout << "Generated type 2 seeds in " << timer.getDuration() << " seconds" << endl;
    }
    
    //Write statistics
    stats.setBinaryFile(extfrag_bin);
    stats.writeStatisticstoFile(outDir, "extended_fragments", num_sampled);
    
    if (hist_path != "") {
        stats.setBinaryFile(type1_bin);
        stats.writeStatisticstoFile(outDir, type1_name, num_sampled);
        
        stats.setBinaryFile(type2_bin);
        stats.writeStatisticstoFile(outDir, type2_name, num_sampled);
    }
    
    string lcloud_out;
    cout << "Writing a line cloud file (TERM Extension) for visualization" << endl;
    secondaryStructureClassifier classifier;
    lcloud_out = outDir + "termext_seed_ca.lcloud";
    fstream out;
    MstUtils::openFile(out, lcloud_out, fstream::out);
    classifier.writeCaInfotoLineFile(extfrag_bin, num_sampled, out);
    out.close();
    
    if (hist_path != "") {
        cout << "Writing a line cloud file (type 1) for visualization" << endl;
        lcloud_out = outDir + "type1_seed_ca.lcloud";
        MstUtils::openFile(out, lcloud_out, fstream::out);
        classifier.writeCaInfotoLineFile(type1_bin, num_sampled, out);
        out.close();
        
        cout << "Writing a line cloud file (type 2) for visualization" << endl;
        lcloud_out = outDir + "type2_seed_ca.lcloud";
        MstUtils::openFile(out, lcloud_out, fstream::out);
        classifier.writeCaInfotoLineFile(type2_bin, num_sampled, out);
        out.close();
    }
    
    
    //Map the seed coverage
    /*
     Vary parameters that filter the seeds
     - match RMSD. The RMSD of the match that generated the seed (depends on the complexity of the fragment)
     - sequence constraint. When applied, only matches with the same residue at the central position are accepted
     */
    
    IC.writePeptideResidues(covDir);
    IC.writeContacts(covDir);
    
    IC.setSeeds(extfrag_bin);
    cout << "Search for segments of seed chains (TERM Extension) that map to the peptide..." << endl;
    IC.findCoveringSeeds();
    cout << "Write coverage to files..." << endl;
    IC.writeAllAlignedSeedsInfo(covDir+"termext_");
    IC.writeBestAlignedSeeds(covDir+"termext_",1,true);
    
    if (hist_path != "") {
        IC.setSeeds(type1_bin);
        cout << "Search for segments of seed chains (type 1) that map to the peptide..." << endl;
        IC.findCoveringSeeds();
        cout << "Write coverage to files..." << endl;
        IC.writeAllAlignedSeedsInfo(covDir+"type1_");
        IC.writeBestAlignedSeeds(covDir+"type1_",1,true);
        
        IC.setSeeds(type2_bin);
        cout << "Search for segments of seed chains (type 2) that map to the peptide..." << endl;
        IC.findCoveringSeeds();
        cout << "Write coverage to files..." << endl;
        IC.writeAllAlignedSeedsInfo(covDir+"type2_");
        IC.writeBestAlignedSeeds(covDir+"type2_",1,true);
    }
    
    // leave commented until future commit where the rest of the pathFromCoveringSeeds is introduced
//    /*
//     Lastly, find the set optimally covering seeds and output as a path string. These can be fused later.
//
//     This is obtained by the following algorithm
//     1) Define the peptide residues to be covered
//     2) Find the seed that covers the most peptide residues. If there is a tie between two seeds, choose
//     the seed with the lowest RMSD. Note: the terminal residues of the seed do not count as covering the peptide
//     unless there is some other seed residue also covering that peptide residue. This is because junctions
//     between seeds *must* have overlap with other seeds (or else they could not be sampled by our method)
//     3) If all peptide residues are covered, terminate. Otherwise, return to step 2
//
//     */
//
//    pathFromCoveringSeeds generatePath(&IC);
//
//    string pathString = generatePath.getCoveringPath();
//
//    if (pathString == "") cout << "Was not able to generate a covering path" << endl;
//    else cout << "path string: " << pathString;
//
//    generatePath.printCoveringSeeds();
    
    cout << "done" << endl;
    return 0;
}
