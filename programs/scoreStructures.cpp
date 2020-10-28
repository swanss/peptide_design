//
//  scoreSeeds.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 4/8/19.
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
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Performs seed scoring using the sequence compatibility score, using seed adjacency graphs for improved performance.");
    opts.addOption("mode", "0 = score seeds, 1 = count contacts for each target residue", false);
    opts.addOption("target", "The target PDB structure", true);
    opts.addOption("structures", "Directory with structures to be scored", true);
    opts.addOption("contacts", "Directory into which to write contact counts, if in contact counting mode, or from which to read the counts if in scoring mode", false);
    opts.addOption("out", "Directory into which to write seed scores", false);
    opts.addOption("config", "The path to a configfile", true);
    opts.addOption("worker", "The worker number (from 1 to numWorkers)", false);
    opts.addOption("numWorkers", "The number of workers", false);
    opts.addOption("numSeedFlank", "Number of residues on either side of the seed residue to use for scoring (default 2)", false);
    opts.addOption("numTargetFlank", "Number of residues on either side of the target residue to use for scoring (default 2)", false);
//    opts.addOption("append", "If true, create new files that cover residues not originally scored by this chunk");
    opts.setOptions(argc, argv);
    
    int mode = opts.getInt("mode", 0);
    string targetPath = opts.getString("target");
//    string graphsPath = opts.getString("graphs");
    string structuresPath = opts.getString("structures");
    string structuresList = MstSystemExtension::join(structuresPath,"structures.list");
    string contactsPath = opts.getString("contacts", "");
    string outPath = opts.getString("out", "");
    
    if (mode == 0 && outPath.size() == 0) MstUtils::error("Need an output path ('--out') for scoring seeds");
    else if (mode == 1 && contactsPath.size() == 0) MstUtils::error("Need a contacts output path ('--contacts') for counting contacts");
    
    string configFile = opts.getString("config");
    int numSeedFlank = opts.getInt("numSeedFlank", 2);
    int numTargetFlank = opts.getInt("numTargetFlank", 2);
    
    int workerIndex = opts.getInt("worker", 1);
    int numWorkers = opts.getInt("numWorkers", 1);
    if (workerIndex < 1 || workerIndex > numWorkers) MstUtils::error("Batch index must be between 1 and numWorkers");
    cout << "Batch " << workerIndex << " of " << numWorkers << endl;
    
    bool append = opts.isGiven("append"); // Only supported for regular scoring, not contact counting
    
    // Initialize
    Structure *target = new Structure(targetPath);
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceCompatibilityScorer scorer(target, rParams, cParams, configFile, numTargetFlank, numSeedFlank, 0.4, 0.05, 0.25, 1, 8000, 0.3);

    if (mode == 1) {
//        // Count contacts only
//        scorer.noQueries = true;
//
//        if (!MstSys::fileExists(contactsPath))
//            MstSys::cmkdir(contactsPath);
//
//        int i = 0;
//        unordered_set<string> scoredSeeds;
//        for (int taskIdx = workerIndex - 1; taskIdx < tasks.size(); taskIdx += numWorkers) {
//            string graphFile = tasks[taskIdx];
//            cout << "graph file: " << graphFile << endl;
//            // Load adjacency graph
//            SeedGraphMap<mstreal> graph;
//            graph.read(MstSystemExtension::join(graphsPath, graphFile), dataPath);
//            scorer.adjacencyGraph = &graph;
//            StructureCache *structures = graph.getStructures();
//
//            for (auto it = structures->begin(); it != structures->end(); it++) {
//                Structure *seed = *it;
//                if (scoredSeeds.count(seed->getName()) != 0) {
//                    cout << "redundant seed!" << endl;
//                    continue;
//                }
//                scoredSeeds.insert(seed->getName());
//
//                if ((i++) % 100 == 0)
//                    cout << "Structure " << i << endl;
//                if (seed->chainSize() > 1)
//                    continue;
//                cout << "Loaded structure " << seed->getName() << endl;
//
//                // Collect contacts from the seed
//                scorer.collectContacts(seed);
//            }
//        }
//
//        cout << "Done counting contacts!" << endl;
//        scorer.writeContactCounts(MstSystemExtension::join(contactsPath, "contact_scores_" + to_string(workerIndex) + ".csv"));
    } else if (mode == 0) {
        // Score seeds
        
        MstSys::cmkdir(outPath);

        // Appending condition
        unordered_set<string> scoredSeeds;
        unordered_map<string, mstreal> preScoredResidues;
//        if (append) {
//            string tempOutPath = MstSystemExtension::join(outPath, "seed_scores_" + to_string(workerIndex) + ".csv");
//            while (MstSys::fileExists(tempOutPath)) {
//                vector<string> lines = MstUtils::fileToArray(tempOutPath);
//                for (string line: lines) {
//                    vector<string> comps = splitString(line, ",");
//                    scoredSeeds.insert(comps[0].substr(0, comps[0].find(":")));
//                    preScoredResidues[comps[0]] = stod(comps[1]);
//                }
//                workerIndex += numWorkers;
//                tempOutPath = MstSystemExtension::join(outPath, "seed_scores_" + to_string(workerIndex) + ".csv");
//            }
//            cout << "New worker index " << workerIndex << ", " << scoredSeeds.size() << " seeds already scored" << endl;
//        }
        
        string fragmentScoresPath = MstSystemExtension::join(outPath, "all_fragment_scores");
        if (!MstSys::fileExists(fragmentScoresPath))
            MstSys::cmkdir(fragmentScoresPath);
        string scoreWritePath = MstSystemExtension::join(fragmentScoresPath, "frag_scores_" + to_string(workerIndex) + ".csv");
        scorer.scoresWritePath = &scoreWritePath;
        
        // Read contact counts
//        if (contactsPath.size() > 0) {
//            int contactFileIdx = 1;
//            string path = MstSystemExtension::join(contactsPath, "contact_scores_" + to_string(contactFileIdx++) + ".txt");
//            while (MstSystemExtension::fileExists(path)) {
//                scorer.readContactCounts(path);
//                path = MstSystemExtension::join(contactsPath, "contact_scores_" + to_string(contactFileIdx++) + ".txt");
//            }
//        }

        double totalTime = 0.0;
        int numTimes = 0;
        
        ofstream outputSS(MstSystemExtension::join(outPath,"structure_scores_") + to_string(workerIndex) + ".tsv", ios::out);
        if (!outputSS.is_open()) {
            cerr << "couldn't open out stream" << endl;
            return 1;
        }
        outputSS << "structure_name\tresidue_name\tscore" << endl;
        
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
            if (scoredSeeds.count(s->getName()) != 0) {
                continue;
            }
            scoredSeeds.insert(s->getName());
            
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
            
            //outputSeedFile.write(seed->getName(), seed->getChain(0).getID(), { writePerResidueScores(result) });
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
        
        cout << "Done! Average time per structures: " << totalTime << "/" << numTimes << " = " << totalTime / numTimes << endl;
        cout << "Scoring statistics: " << endl;
        cout << "\t" << scorer.numResiduesScored() << " residues scored" << endl;
        cout << "\t" << scorer.numUniqueResiduesScored() << " unique residues scored" << endl;
        cout << "\t" << scorer.numSeedQueries() << " seed queries" << endl;
        cout << "\t" << scorer.numTargetQueries() << " target queries" << endl;

    }
    
    delete target;

    return 0;
}
