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
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Scores the interface formed between a set of peptides and the target protein.");
    opts.addOption("target", "The target PDB structure", false);
    opts.addOption("structures", "Directory with structures to be scored", false);
    opts.addOption("list", "File with names of structures to be scored."); 
    opts.addOption("contacts", "Directory into which to write contact counts, if in contact counting mode, or from which to read the counts if in scoring mode", false);
    opts.addOption("out", "Directory into which to write seed scores", false);
    opts.addOption("config", "The path to a configfile", true);
    opts.addOption("worker", "The worker number (from 1 to numWorkers)", false);
    opts.addOption("numWorkers", "The number of workers", false);
    opts.addOption("numSeedFlank", "Number of residues on either side of the seed residue to use for scoring (default 2)", false);
    opts.addOption("numTargetFlank", "Number of residues on either side of the target residue to use for scoring (default 2)", false);
    opts.addOption("complex", "The path to a peptide-protein complex PDB (single complex mode)",false);
    opts.addOption("peptide", "The peptide chain name (single complex mode)",false);
    opts.addOption("dtermenConfig", "Path to dTERMen config file.");
    opts.addOption("filter", "Inidcator of whether to filter on freedom.", false);
    opts.addOption("freedomLim", "Upper bound on freedom to use if filtering.");
    opts.addOption("scoreType", "Indicates which scoring function to use.");
    opts.setOptions(argc, argv);
    
    if ((!opts.isGiven("target")|!opts.isGiven("structures"))&(!opts.isGiven("complex")|!(opts.isGiven("peptide")))) MstUtils::error("Must provide either --target/--structures or --complex/--peptide","scoreStructures");
    
    int mode = opts.getInt("mode", 0);
    string targetPath = opts.getString("target");
    string structuresPath = opts.getString("structures");
    string structuresList;
    if ((!opts.isGiven("list")))
        structuresList = MstSystemExtension::join(structuresPath,"structures.list");
    else
        structuresList = opts.getString("list");
    string contactsPath = opts.getString("contacts", "");
    string outPath = opts.getString("out", "");
    
    string complexPath = opts.getString("complex","");
    string peptideChainID = opts.getString("peptide","");

    string dtermenConfig = opts.getString("dtermenConfig");
    bool filter = opts.getBool("filter");
    mstreal freedomLim = opts.getReal("freedomLim");
    int scoreType = opts.getInt("scoreType");
    if ((scoreType < 1) || (scoreType > 4)) MstUtils::error("scoreType must be 1-4.","scoreStructures");

     
    if (mode == 0 && outPath.size() == 0) MstUtils::error("Need an output path ('--out') for scoring seeds");
    else if (mode == 1 && contactsPath.size() == 0) MstUtils::error("Need a contacts output path ('--contacts') for counting contacts");

    string subsetOutputPath = opts.getString("subsetOutputPath");
    string peptideStructurePath = opts.getString("peptideStructurePath");
    
    string configFile = opts.getString("config");
    int numSeedFlank = opts.getInt("numSeedFlank", 2);
    int numTargetFlank = opts.getInt("numTargetFlank", 2);
    
    int workerIndex = opts.getInt("worker", 1);
    int numWorkers = opts.getInt("numWorkers", 1);
    if (workerIndex < 1 || workerIndex > numWorkers) MstUtils::error("Batch index must be between 1 and numWorkers");
    cout << "Batch " << workerIndex << " of " << numWorkers << endl;
        
    // Initialize
    Structure *target = nullptr;
    Structure *peptide = nullptr;
    if (complexPath != "") {
        target = new Structure(complexPath);
        Chain* pChain = target->getChainByID(peptideChainID);
        peptide = new Structure(*pChain);
        peptide->setName(MstSystemExtension::fileName(complexPath));
        target->deleteChain(pChain);
    } else {
        target = new Structure(targetPath);
    }
    
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceStructureCompatibilityScorer scorer(target, rParams, cParams, configFile, dtermenConfig, filter, freedomLim, scoreType, numTargetFlank, numSeedFlank, 0.4, 0.005, 0.25, 1, 8000, 0.3);

    // Score seeds
    MstSys::cmkdir(outPath);
    
    string fragmentScoresPath = MstSystemExtension::join(outPath, "all_fragment_scores");
    if (!MstSys::fileExists(fragmentScoresPath))
        MstSys::cmkdir(fragmentScoresPath);
    string scoreWritePath = MstSystemExtension::join(fragmentScoresPath, "frag_scores_" + to_string(workerIndex) + ".csv");
    string filterWritePath = MstSystemExtension::join(fragmentScoresPath, "filtered_res_" + to_string(workerIndex) + ".csv");
    string freedomWritePath = MstSystemExtension::join(fragmentScoresPath, "freedom_info_" + to_string(workerIndex) + ".csv");
    scorer.scoresWritePath = &scoreWritePath;
    scorer.filterWritePath = &filterWritePath;
    scorer.freedomWritePath = &freedomWritePath;
    scorer.maxScore = 100;

    
    double totalTime = 0.0;
    int numTimes = 0;
    
    ofstream outputSS(MstSystemExtension::join(outPath,"structure_scores_") + to_string(workerIndex) + ".tsv", ios::out);
    if (!outputSS.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return 1;
    }
    outputSS << "structure_name\tresidue_name\tscore" << endl;
    
    if (peptide != nullptr) {
        cout << "Loaded structure " << peptide->getName() << endl;
        
        // Score the seed
        high_resolution_clock::time_point startTime = high_resolution_clock::now();
        auto result = scorer.score(peptide);
        high_resolution_clock::time_point endTime = high_resolution_clock::now();
        double time = duration_cast<seconds>( endTime - startTime ).count();
        cout << "time: " << time << endl;
        totalTime += time;
        numTimes += 1;
        
        //outputSeedFile.write(seed->getName(), seed->getChain(0).getID(), { writePerResidueScores(result) });
        for (auto resScore: result) {
            outputSS << peptide->getName() << "\t" << resScore.first->getChainID() << resScore.first->getNum() << "\t" << resScore.second << endl;
        }
    } else  {
        cout << "loading structure cache. List: " << structuresList << " Prefix: " << structuresPath << endl;
        StructureCache* structures = new StructureCache(structuresPath);
        structures->preloadFromPDBList(structuresList);
        long cache_size = structures->size();
        
        int i = 0;
        auto it = structures->begin();
        for (int i = 0; i < workerIndex-1; i++) {
            it++;
            if (it == structures->end()) break;
        }
        while (it != structures->end()) {
            Structure *s = *it;
            
            if ((i++) % 100 == 0)
                cout << "Structure " << i << endl;
            if (s->chainSize() > 1)
                continue;
            cout << "Loaded structure " << s->getName() << endl;
            
            // Score the seed
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            auto result = scorer.score(s);
            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            double time = duration_cast<seconds>( endTime - startTime ).count();
            cout << "time: " << time << endl;
            totalTime += time;
            numTimes += 1;
            
            for (auto resScore: result) {
                outputSS << s->getName() << "\t" << resScore.first->getChainID() << resScore.first->getNum() << "\t" << resScore.second << endl;
            }
            
            // advance the iterator
            for (int i = 0; i < numWorkers; i++) {
                it++;
                if (it == structures->end()) break;
            }
        }
        delete structures;
    }
    
    cout << "Done! Average time per structures: " << totalTime << "/" << numTimes << " = " << totalTime / numTimes << endl;
    cout << "Scoring statistics: " << endl;
    cout << "\t" << scorer.numResiduesScored() << " residues scored" << endl;
    cout << "\t" << scorer.numUniqueResiduesScored() << " unique residues scored" << endl;
    cout << "\t" << scorer.numSeedQueries() << " seed queries" << endl;
    cout << "\t" << scorer.numTargetQueries() << " target queries" << endl;

    if (peptide != nullptr) delete peptide;
    delete target;

    return 0;
}
