#include <cstdio>

//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "utilities.h"
#include "coverage.h"
#include "benchmarkutilities.h"
//#include "termextension.h"
#include "secondarystructure.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Maps seeds to the peptide in a peptide-protein complex and computes various statistics. At the end, a minimum covering set is determined and these seeds are fused together.");
    op.addOption("bin_path", "path to a seed binary file.",true);
    op.addOption("pdb", "Structure file containing peptide chain and protein chain(s)", true);
    op.addOption("peptide", "Peptide chain ID", true);
    op.addOption("max_rmsd", "The max RMSD threshold used when determining whether a seed aligns to the peptide or not. (default 2.0)");
    op.addOption("max_match_num", "The max match number (matches are ranked by RMSD) that is accepted when considering seeds. (default 100000)");
    op.addOption("config","Path to the configuration file (specifies fasst database and rotamer library)",true);
    op.addOption("hist","Path to the histogram file with the seed distance distribution that will be matched in the null model seeds");
    op.addOption("write_all_aligned","If provided, write all seeds that align to the peptide with rmsd < max_rmsd. Omitting this option speeds up the runtime and saves space");
    op.addOption("max_seed_length_fuse","The maximum length of seed segment we should try to use for fusing. The algorithm will try to search for seeds of this length first and only try using shorter seed segments if there is nothing at this length with with rmsd < max_rmsd. The minimum length of seed segments is 3. (default 5)");
    op.addOption("force_chimera","If provided, then the seed segments that are selected as covering the peptide must all be sourced from different seeds.");
    op.addOption("two_step_fuse","When fusing, optimize without and then with internal coordinate constraints.");
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
    int max_match_num = op.getReal("max_match_num",100000);
    configFile config(op.getString("config"));
    string hist_path = op.getString("hist","");
    bool write_all_aligned = op.isGiven("write_all_aligned");
    int max_seed_length_fuse = op.getInt("max_seed_length_fuse",5);
    bool force_chimera = op.isGiven("force_chimera");
    bool two_step_fuse = op.isGiven("two_step_fuse");
    
    string complex_name = MstSys::splitPath(complex.getName(), 1);
    string pdb_id = complex_name.substr(0,4); //PDB ID is always 4 characters
    
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
    string fusDir = "fuse_coverage_output/";
    MstSys::cmkdir(outDir,makeParents);
    MstSys::cmkdir(covDir,makeParents);
    MstSys::cmkdir(fusDir,makeParents);
    
    //Initialize the interface coverage class
    interfaceCoverage IC(complex, p_cid, config.getRL());
    IC.setMaxRMSD(max_rmsd);
    IC.setMaxSegmentLength(max_seed_length);
    IC.setMaxMatchNumber(max_match_num);
    
    // remove the chain from a copy of the complex
    Structure target(complex);
    Chain* peptideChain = target.getChainByID(p_cid);
    target.deleteChain(peptideChain);
        
    //Randomize seeds
    /**
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
        oneDimBinnedData proposal = stats.generateDistanceHistogram();
        string proposal_hist_path = outDir + "seed_centroid_distance_proposal.csv";
        proposal.writeHistFile(proposal_hist_path);
        
        //generate type 2 null model seeds with rejection sampling
        oneDimBinnedData target(hist_path);
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
    
    
    //Map the seed coverage and, if possible, fuse covering seeds
    /**
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
    if (write_all_aligned) IC.writeAllAlignedSeedsInfo(covDir+"termext_");
    IC.writeBestAlignedSeeds(covDir+"termext_",1,true);
    
    coverageBenchmarkUtils::fuseCoveringSeeds(&IC, force_chimera, max_seed_length_fuse, fusDir+"termext", pdb_id, target, complex, two_step_fuse, config.getRL());
    
    if (hist_path != "") {
        IC.setSeeds(type1_bin);
        cout << "Search for segments of seed chains (type 1) that map to the peptide..." << endl;
        IC.findCoveringSeeds();
        cout << "Write coverage to files..." << endl;
        if (write_all_aligned) IC.writeAllAlignedSeedsInfo(covDir+"type1_");
        IC.writeBestAlignedSeeds(covDir+"type1_",1,true);

        coverageBenchmarkUtils::fuseCoveringSeeds(&IC, force_chimera, max_seed_length_fuse, fusDir+"type1", pdb_id, target, complex, two_step_fuse, config.getRL());

        
        IC.setSeeds(type2_bin);
        cout << "Search for segments of seed chains (type 2) that map to the peptide..." << endl;
        IC.findCoveringSeeds();
        cout << "Write coverage to files..." << endl;
        if (write_all_aligned) IC.writeAllAlignedSeedsInfo(covDir+"type2_");
        IC.writeBestAlignedSeeds(covDir+"type2_",1,true);
        
        coverageBenchmarkUtils::fuseCoveringSeeds(&IC, force_chimera, max_seed_length_fuse, fusDir+"type2", pdb_id, target, complex, two_step_fuse, config.getRL());
    }
    
    // Delete the random baseline seed binary files (they could be regenerated later)
    char *type1_bin_char = &type1_bin[0];
    remove(type1_bin_char);
    
    char *type2_proposal_bin_char = &type2_proposal_bin[0];
    remove(type2_proposal_bin_char);
    
    char *type2_bin_char = &type2_bin[0];
    remove(type2_bin_char);
    
    cout << "Done" << endl;
    return 0;
}
