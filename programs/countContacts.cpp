//
//  countContacts.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 1/6/21.
//  Copyright © 2021 Sebastian Swanson. All rights reserved.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "utilities.h"
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    /*
     This program counts contacts between seeds and the target. This is to calculate the normalizing
     constant that is used later to find the "seed score", which is designed to reflect the probability of
     the seed, given the target. This is calculated for every contact between a seed and target residue.
     
     Seed score -> P(seed residue | target residue) = # of overlapping seeds / # of seeds that contact
     the target residue
     
     Since the maximum number of possible overlapping seeds is the sum of all length-k seed windows,
     we consider all contacts between target residues and seed residues *except* seed residues within
     floor(k/2) from the terminus. These "terminal" residues do not have sufficient flanking residues
     to be considered for overlaps and thus should not be considered in defining the normalizing constant
     
     All types of contacts are considered (sidechain-sidechain, sidechain-backbone, backbone-sidechain,
     and backbone-backbone), but each target residue/seed residue pair can be defined as in contact
     no more than once.
     */
 
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Counts the number of contacts between protein residues and seeds and reports this as contacts per protein residue.");
    opts.addOption("target", "The target PDB structure", false);
    opts.addOption("base","The name that shold be appended to the output filename (inferred in concat mode)",false);
    opts.addOption("seed_out","If provided, will also output a .txt file with the number of contacts per each seed residue");
    opts.addOption("bin", "A seed binary file containing all seeds", false);
    opts.addOption("configFile","A configuration file specifying the rotamer library/fasstDB",false);
    opts.addOption("flankRes","If provided, only counts contacts that have this many flanking residues on each side",false);
    opts.addOption("worker", "The worker number (from 1 to numWorkers)", false);
    opts.addOption("numWorkers", "The number of workers", false);
    opts.addOption("combine","A file where each line is a contacts file. Contact counts from each file are combined into a single file.",false);
    opts.setOptions(argc, argv);
    
    if ((!opts.isGiven("combine")) && ((!opts.isGiven("target")) || (!opts.isGiven("base")) || (!opts.isGiven("bin")) || (!opts.isGiven("configFile")))) MstUtils::error("If not in combine mode, then --target, --base, --bin and --combine must be provided");
    
    string targetPath = opts.getString("target");
    string baseName = opts.getString("base");
    string seedBin = opts.getString("bin");
    string configFile = opts.getString("configFile");
    int flankRes = opts.getInt("flankRes",0);
    int workerIndex = opts.getInt("worker",1);
    int numWorkers = opts.getInt("numWorkers",1);
    
    class configFile configObj(configFile);
    
    RotamerLibrary RL;
    RL.readRotamerLibrary(configObj.getRL());
    
    Structure *target = new Structure(targetPath);
    contactParams cParams(3.5,0.01,0.01);
    
    contactCounter cCounter(target,&RL,cParams,flankRes);
    
    if (opts.isGiven("seed_out")) {
        string allSeedContactsFile = baseName + "_" + MstUtils::toString(workerIndex) + "_" + "seedConts.tsv";
        cCounter.setWriteSeedContacts(allSeedContactsFile);
    }
    
    if (opts.isGiven("combine")) {
        
        vector<string> contactsFilesList = MstUtils::fileToArray(opts.getString("combine"));
        
        for (string contactsFile : contactsFilesList) {
            cCounter.readContactsFile(contactsFile);
        }
        
        string contactsFile = baseName + "_" + "conts.tsv";
        cCounter.writeContactsFile(contactsFile);
        
    } else {
        int batchSize = 1000;
        string chainID = "0";
        StructureIterator structIter = StructureIterator(seedBin,batchSize,chainID,workerIndex,numWorkers);
        
        while (structIter.hasNext()) {
            vector<Structure*> seed_chains = structIter.next();
            for (Structure* s : seed_chains) {
                // count the contacts between this seed structure and the target protein
                cCounter.countContacts(s);
            }
        }
        string contactsFile = baseName + "_" + MstUtils::toString(workerIndex) + "_" + "conts.tsv";
        cCounter.writeContactsFile(contactsFile);
    }
    
    delete target;
    
    return 0;
}
