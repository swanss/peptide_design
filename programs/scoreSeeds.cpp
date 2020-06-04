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
    opts.addOption("graphs", "Directory in which seed graphs are stored", true);
    opts.addOption("data", "Directory in which seed structures are stored", true);
    opts.addOption("contacts", "Directory into which to write contact counts, if in contact counting mode, or from which to read the counts if in scoring mode", false);
    opts.addOption("out", "Directory into which to write seed scores", false);
    opts.addOption("fasstDB", "The target PDB files to be used in scoring searches.", false);
    opts.addOption("batch", "The batch number (from 1 to numBatches)", false);
    opts.addOption("numBatches", "The number of batches", false);
    opts.addOption("numSeedFlank", "Number of residues on either side of the seed residue to use for scoring (default 2)", false);
    opts.addOption("numTargetFlank", "Number of residues on either side of the target residue to use for scoring (default 2)", false);
    opts.addOption("append", "If true, create new files that cover residues not originally scored by this chunk");
    opts.setOptions(argc, argv);
    
    int mode = opts.getInt("mode", 0);
    string targetPath = opts.getString("target");
    string graphsPath = opts.getString("graphs");
    string dataPath = opts.getString("data");
    string contactsPath = opts.getString("contacts", "");
    string outPath = opts.getString("out", "");
    if (mode == 0 && outPath.size() == 0) {
        cerr << "Need an output path for scoring seeds" << endl;
        return 1;
    } else if (mode == 1 && contactsPath.size() == 0) {
        cerr << "Need a contacts output path for counting contacts" << endl;
        return 1;
    }
    
    string fasstDB = opts.getString("fasstDB", "/home/grigoryanlab/library/databases/dTERMen.databases/2019-01-22/dtermen.sim");
    int numSeedFlank = opts.getInt("numSeedFlank", 2);
    int numTargetFlank = opts.getInt("numTargetFlank", 2);
    
    int batchIndex = opts.getInt("batch", 1);
    int numBatches = opts.getInt("numBatches", 1);
    if (batchIndex < 1 || batchIndex > numBatches) {
        cerr << "Batch index must be between 1 and numBatches" << endl;
        return 1;
    }
    cout << "Batch " << batchIndex << " of " << numBatches << endl;
    bool append = opts.isGiven("append"); // Only supported for regular scoring, not contact counting
    
    // Initialize
    Structure *target = new Structure(targetPath);
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceCompatibilityScorer scorer(target, rParams, cParams, fasstDB, numTargetFlank, numSeedFlank, 0.4, 0.05, 0.25, 1, 8000, 0.3);

    vector<string> tasks = MstUtils::fileToArray(graphsPath + "tasks.txt");

    if (mode == 1) {
        // Count contacts only
        scorer.noQueries = true;
        
        if (!MstSys::fileExists(contactsPath))
            MstSys::cmkdir(contactsPath);
        
        int i = 0;
        unordered_set<string> scoredSeeds;
        for (int taskIdx = batchIndex - 1; taskIdx < tasks.size(); taskIdx += numBatches) {
            string graphFile = tasks[taskIdx];
            cout << "graph file: " << graphFile << endl;
            // Load adjacency graph
            SeedGraphMap<mstreal> graph;
            graph.read(MstSystemExtension::join(graphsPath, graphFile), dataPath);
            scorer.adjacencyGraph = &graph;
            StructureCache *structures = graph.getStructures();
            
            for (auto it = structures->begin(); it != structures->end(); it++) {
                Structure *seed = *it;
                if (scoredSeeds.count(seed->getName()) != 0) {
                    cout << "redundant seed!" << endl;
                    continue;
                }
                scoredSeeds.insert(seed->getName());
                
                if ((i++) % 100 == 0)
                    cout << "Structure " << i << endl;
                if (seed->chainSize() > 1)
                    continue;
                cout << "Loaded structure " << seed->getName() << endl;
                
                // Collect contacts from the seed
                scorer.collectContacts(seed);
            }
        }
        
        cout << "Done counting contacts!" << endl;
        scorer.writeContactCounts(MstSystemExtension::join(contactsPath, "contact_scores_" + to_string(batchIndex) + ".csv"));
    } else if (mode == 0) {
        // Score seeds
        
        if (!MstSys::fileExists(outPath))
            MstSys::cmkdir(outPath);

        // Appending condition
        unordered_set<string> scoredSeeds;
        unordered_map<string, mstreal> preScoredResidues;
        if (append) {
            string tempOutPath = MstSystemExtension::join(outPath, "seed_scores_" + to_string(batchIndex) + ".csv");
            while (MstSys::fileExists(tempOutPath)) {
                vector<string> lines = MstUtils::fileToArray(tempOutPath);
                for (string line: lines) {
                    vector<string> comps = splitString(line, ",");
                    scoredSeeds.insert(comps[0].substr(0, comps[0].find(":")));
                    preScoredResidues[comps[0]] = stod(comps[1]);
                }
                batchIndex += numBatches;
                tempOutPath = MstSystemExtension::join(outPath, "seed_scores_" + to_string(batchIndex) + ".csv");
            }
            cout << "New batch index " << batchIndex << ", " << scoredSeeds.size() << " seeds already scored" << endl;
        }
        
        string fragmentScoresPath = MstSystemExtension::join(outPath, "all_fragment_scores");
        if (!MstSys::fileExists(fragmentScoresPath))
            MstSys::cmkdir(fragmentScoresPath);
        string scoreWritePath = MstSystemExtension::join(fragmentScoresPath, "frag_scores_" + to_string(batchIndex) + ".csv");
        scorer.scoresWritePath = &scoreWritePath;
        
        // Read contact counts
        if (contactsPath.size() > 0) {
            int contactFileIdx = 1;
            string path = MstSystemExtension::join(contactsPath, "contact_scores_" + to_string(contactFileIdx++) + ".txt");
            while (MstSystemExtension::fileExists(path)) {
                scorer.readContactCounts(path);
                path = MstSystemExtension::join(contactsPath, "contact_scores_" + to_string(contactFileIdx++) + ".txt");
            }
        }

        double totalTime = 0.0;
        int numTimes = 0;
        
        ofstream outputSS(outPath + "seed_scores_" + to_string(batchIndex) + ".csv", ios::out);
        if (!outputSS.is_open()) {
            cerr << "couldn't open out stream" << endl;
            return 1;
        }
        
        int i = 0;
        for (int taskIdx = batchIndex - 1; taskIdx < tasks.size(); taskIdx += numBatches) {
            string graphFile = tasks[taskIdx];
            cout << "graph file: " << graphFile << endl;
            // Load adjacency graph
            SeedGraphMap<mstreal> graph;
            graph.read(MstSystemExtension::join(graphsPath, graphFile), dataPath);
            
            // Add in values from pre-scored residues
            for (auto item: preScoredResidues) {
                Residue *res = graph.getResidueFromFile(item.first, false);
                if (res != nullptr && graph.contains(res)) {
                    graph.setValue(res, item.second);
                }
            }
            
            scorer.adjacencyGraph = &graph;
            StructureCache *structures = graph.getStructures();
            
            for (auto it = structures->begin(); it != structures->end(); it++) {
                Structure *seed = *it;
                if (scoredSeeds.count(seed->getName()) != 0) {
                    continue;
                }
                scoredSeeds.insert(seed->getName());
                
                if ((i++) % 100 == 0)
                    cout << "Structure " << i << endl;
                if (seed->chainSize() > 1)
                    continue;
                cout << "Loaded structure " << seed->getName() << endl;
                
                // Score the seed
                high_resolution_clock::time_point startTime = high_resolution_clock::now();
                auto result = scorer.score(seed);
                high_resolution_clock::time_point endTime = high_resolution_clock::now();
                double time = duration_cast<seconds>( endTime - startTime ).count();
                cout << "time: " << time << endl;
                totalTime += time;
                numTimes += 1;
                
                //outputSeedFile.write(seed->getName(), seed->getChain(0).getID(), { writePerResidueScores(result) });
                for (auto resScore: result) {
                    outputSS << graph.writeCodeForResidue(resScore.first) << "," << resScore.second << endl;
                }
            }
        }
        
        cout << "Done! Average time per seed: " << totalTime << "/" << numTimes << " = " << totalTime / numTimes << endl;
        cout << "Scoring statistics: " << endl;
        cout << "\t" << scorer.numResiduesScored() << " residues scored" << endl;
        cout << "\t" << scorer.numUniqueResiduesScored() << " unique residues scored" << endl;
        cout << "\t" << scorer.numSeedQueries() << " seed queries" << endl;
        cout << "\t" << scorer.numTargetQueries() << " target queries" << endl;

    }
    
    delete target;

    return 0;
}
