#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "fusehelper.h"
#include "seedgraph.h"
#include "seedscore.h"
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    string basePath = "/home/ifs-users/venkats/dtermen/binding_site_new_seeds/";
    string dataPath = "/home/ifsdata/scratch/grigoryanlab/swans/peptide_design_benchmark/seed_retrieval/target_centered_approach/retrieveSeedsTest3/1DUY/output/seeds/";
    string outPath = basePath + "seed_scores_native/";
    string fasstDB = "/home/grigoryanlab/library/databases/dTERMen.databases/2019-01-22/dtermen.sim";

    Structure *target = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb");
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceCompatibilityScorer scorer(target, rParams, cParams, fasstDB, 2, 2, 0.4, 0.05, 0.25, 1, 8000, 0.7); // new StructureCompatibilityScorer(target, fParams, rParams, cParams, fasstDB);
    string scoreWritePath = outPath + "frag_scores.csv";
    scorer.scoresWritePath = &scoreWritePath;

    // Read contact counts
    string contactCountsPath = basePath + "seed_scores_no_clash/contact_scores/";
    for (int i = 1; i <= 35; i++) {
        scorer.readContactCounts(contactCountsPath + "contact_scores_" + to_string(i) + ".csv");
    }
    
    ofstream outputSS(outPath + "residue_scores.csv", ios::out);
    if (!outputSS.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return 1;
    }
    
    // Load adjacency graph
    string graphsPath = basePath + "graphs/";
    SeedGraphMap<mstreal> graph;
    graph.read(MstSystemExtension::join(graphsPath, "whole_graph.txt"), dataPath);
    scorer.adjacencyGraph = &graph;
    
    vector<string> seedPaths = MstUtils::fileToArray(outPath + "poses/pdbs.txt");
    for (string seedName: seedPaths) {
        cout << "scoring " << seedName << endl;
        Structure *seed = new Structure(outPath + "poses/" + seedName);
        
        high_resolution_clock::time_point startTime = high_resolution_clock::now();
        auto result = scorer.score(seed);
        high_resolution_clock::time_point endTime = high_resolution_clock::now();
        cout << "Scoring took " << duration_cast<seconds>( endTime - startTime ).count() << " seconds" << endl;
        cout << "Result: ";
        for (auto it: result) {
            cout << it.first->getName() << it.first->getResidueIndex() << " = " << it.second << ", " << endl;
        }
        cout << endl;
        for (auto resScore: result) {
            outputSS << seed->getName() << ":" << resScore.first->getResidueIndex() << "," << resScore.second << endl;
        }
        delete seed;
    }

    delete target;
    
    return 0;
}
