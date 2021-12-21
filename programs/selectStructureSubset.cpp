#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "utilities.h"
#include "makeclustertestset.h"
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Scores the interface formed between a set of peptides and the target protein.");
    opts.addOption("target", "The target PDB structure", false);
    opts.addOption("structures", "Directory with structures to be scored", false);
    opts.addOption("config", "The path to a configfile", true);
    opts.addOption("out", "Indicates name of file to write list of subset of structures to score.", false);
    opts.addOption("peptide", "Path to structure of peptide.");
    opts.addOption("negType", "Indicates how to choose negative controls.");
    opts.addOption("num", "Number of structures to choose.");
    opts.setOptions(argc, argv);
      
    string targetPath = opts.getString("target");
    string structuresPath = opts.getString("structures");
    string structuresList = MstSystemExtension::join(structuresPath,"structures.list");
    
    string subsetOutputPath = opts.getString("out");
    string peptideStructurePath = opts.getString("peptide");
    
    string configFile = opts.getString("config");
    int negType = opts.getInt("negType");

    int nStructures = opts.getInt("num");

    cout << "Neg type: " << negType << endl;
         
    // Initialize
    Structure *target = nullptr;
    Structure *peptide = nullptr;
    target = new Structure(targetPath);
    peptide = new Structure(peptideStructurePath);

    clusterTestSet clusterer(structuresPath, structuresList, target, peptide, configFile, nStructures);
    clusterer.selectStructures(negType);
    cout << "Created subset of structures. Selected " << clusterer.getNumTruePositives() << " true positives, " << clusterer.getNumTrueNegatives() << " true negatives, and " << clusterer.getNumUnkownSamples() << " random structures." << endl;
    cout << "Dist Limit: " << clusterer.getDistLimit() << endl;
    clusterer.writeStructures(subsetOutputPath); 
} 