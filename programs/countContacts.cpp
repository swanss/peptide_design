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
    opts.addOption("structures_list","A file where each line is the path to a PDB file describing a seed",false);
    opts.addOption("bin", "A seed binary file containing all seeds", false);
    opts.addOption("configFile","A configuration file specifying the rotamer library/fasstDB",false);
    opts.addOption("flankRes","If provided, only counts contacts that have this many flanking residues on each side",false);
    opts.addOption("worker", "The worker number (from 1 to numWorkers)", false);
    opts.addOption("numWorkers", "The number of workers", false);
    opts.addOption("combine","A file where each line is a contacts file. Contact counts from each file are combined into a single file.",false);
    
    //The options are used when scoring a single complex
    opts.addOption("complex", "The path to a PDB structure containing the peptide and protein chains (single complex mode)",false);
    opts.addOption("peptide", "The peptide chain ID (single complex mode)",false);
    opts.setOptions(argc, argv);
    
    if (opts.isGiven("complex") and opts.isGiven("peptide")) {
        cout << "PDB of complex and chain name provided, will score single complex" << endl;
    } else if ((!opts.isGiven("combine")) && ((!opts.isGiven("target")) || (!opts.isGiven("base")) || !(opts.isGiven("bin") || opts.isGiven("structures_list")) || (!opts.isGiven("configFile")))) MstUtils::error("If not in combine mode, then --target, --base, and either --bin or --structures_list must be provided");
    
    string targetPath = opts.getString("target");
    string baseName = opts.getString("base");
    string structureList = opts.getString("structures_list","");
    string seedBin = opts.getString("bin","");
    string configFile = opts.getString("configFile");
    int flankRes = opts.getInt("flankRes",0);
    int workerIndex = opts.getInt("worker",1);
    int numWorkers = opts.getInt("numWorkers",1);
    
    string complexPath = opts.getString("complex","");
    string peptideChainID = opts.getString("peptide","");
    
    class configFile configObj(configFile);
    
    RotamerLibrary RL;
    RL.readRotamerLibrary(configObj.getRL());
    
    Structure* peptide = nullptr; //variable only used in single complex mode
    
    Structure *target = nullptr;
    if (complexPath != "") {
        target = new Structure(complexPath);
        Chain* pChain = target->getChainByID(peptideChainID);
        peptide = new Structure(*pChain);
        peptide->setName(MstSystemExtension::fileName(complexPath));
        target->deleteChain(pChain);
    } else target = new Structure(targetPath);
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
        
    } else if (peptide != nullptr) {
        // single complex mode
        cCounter.countContacts(peptide);
        string contactsFile = baseName + "_" + MstUtils::toString(workerIndex) + "_" + "conts.tsv";
        cCounter.writeContactsFile(contactsFile);
    } else {
        // batch scoring mode
        int batchSize = 1000;
        string chainID = "0";
        StructureIterator* structIter = nullptr;
        if (structureList != "") {
            vector<string> structures_paths = MstUtils::fileToArray(structureList);
            cout << "Structure list with " << structures_paths.size() << " structures" << endl;
            vector<string>* chain_ids = nullptr;
            structIter = new StructureIterator(structures_paths,batchSize,chain_ids,workerIndex,numWorkers);
        } else if (seedBin != "") {
            structIter = new StructureIterator(seedBin,batchSize,chainID,workerIndex,numWorkers);
        } else {
            MstUtils::error("Must provide either structure list or seed binary file","countContacts::main");
        }
        
        while (structIter->hasNext()) {
            vector<Structure*> seed_chains = structIter->next();
            for (Structure* s : seed_chains) {
                // count the contacts between this seed structure and the target protein
                cCounter.countContacts(s);
            }
        }
        string contactsFile = baseName + "_" + MstUtils::toString(workerIndex) + "_" + "conts.tsv";
        cCounter.writeContactsFile(contactsFile);
    }
    
    if (peptide == nullptr) delete peptide;
    delete target;
    
    return 0;
}
