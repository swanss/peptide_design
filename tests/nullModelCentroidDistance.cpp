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
    op.addOption("config", "Path to the configuration file (specifies fasst database and rotamer library)",true);
    op.addOption("hist", "Path to the histogram file with the seed distance distribution that will be matched in the null model seeds",true);
    op.addOption("no_clash_check", "Turns off filtering seed poses with clashes to the protein");
    op.addOption("no_rejection_sampling", "Turns off rejection sampling to approximate a target distribution");
    op.setOptions(argc, argv);
    
    MstTimer timer;
    
    // Variables set at the time of compilation
    int num_sampled = INT_MAX; //this is approximately the number of seeds that will be sampled when writing line clouds
    
    // Variables provided by user
    string extfrag_bin = op.getString("bin_path");
    Structure complex(op.getString("pdb"));
    string p_cid = op.getString("peptide");
    configFile config(op.getString("config"));
    string hist_path = op.getString("hist");
    
    // Make directories
    bool makeParents = true;
    // Fragment output folder
    string outDir = "output/";
    MstSys::cmkdir(outDir,makeParents);
    
    //Randomize seeds
    /*
     Two types of randomized seeds
     Type 1) randomized orientation, position sampled from previous seed
     Type 2) randomized position and orientation
     */
    string type2_name = "type2_seeds";
    string type2_bin = outDir + "type2_seeds.bin";
    bool position = true;
    bool orientation = true;
    
    cout << "Generating null model seeds..." << endl;
    
    naiveSeedsFromDB naiveSeeds(complex, p_cid, extfrag_bin, config.getDB(), config.getRL());
    if (op.isGiven("no_clash_check")) naiveSeeds.setClashChecking(false);
    if (op.isGiven("no_rejection_sampling")) naiveSeeds.setRejectionSampling(false);
    
    timer.start();
    naiveSeeds.newPose(outDir, type2_name, position, orientation);
    timer.stop();
    cout << "Generated type 2 seeds in " << timer.getDuration() << " seconds" << endl;
    
    
    //Write statistics
    seedStatistics stats(complex, p_cid, extfrag_bin); //construct with term seeds
    stats.writeStatisticstoFile(outDir, "extended_fragments", num_sampled);
    
    stats.setBinaryFile(type2_bin); //set to type 2 seeds
    stats.writeStatisticstoFile(outDir, type2_name, num_sampled);
    
    cout << "done" << endl;
    return 0;
}
