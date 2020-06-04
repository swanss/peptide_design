//
//  testPathFinder.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 4/22/19.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "findpaths.h"

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    RankList<string> testRanks(3);
    testRanks.insert("a", 5.0);
    testRanks.insert("b", 3.0);
    testRanks.insert("c", 2.0);
    testRanks.insert("d", 8.0);
    vector<string> expectedResult = { "d", "a", "b" };
    if (testRanks.items() != expectedResult) {
        cout << "incorrect rank" << endl;
        for (string val: testRanks.items()) {
            cout << val << ",";
        }
        cout << endl;
    }
    
    string basePath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/1DUY_overlaps_dot_75/";
    string graphFile = basePath + "bond_graphs/whole_graph.txt"; //"/home/ifs-users/venkats/dtermen/binding_site_new_seeds/bond_seed_graphs/whole_graph.txt";
    string seedPath = "/home/ifsdata/scratch/grigoryanlab/swans/term_based_design_method/old_seed_dirs/1DUY/output/seeds/";
    string scoresPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/1DUY_flank2/seed_scores/";

    SeedGraph graph;
    graph.read(graphFile, seedPath);
    cout << graph.residueSize() << " residues" << endl;
    cout << graph.edgeSize() << " edges" << endl;

    SeedGraphMap<mstreal> scoreMap(graph.getStructures());
    scoreMap.read(graphFile, seedPath, &graph);
    cout << "Finished reading map" << endl;
    
    for (int i = 1; i <= 40; i++) {
        vector<string> lines = MstUtils::fileToArray(scoresPath + "seed_scores_" + to_string(i) + ".csv");
        for (string line: lines) {
            vector<string> comps = splitString(line, ",");
            Residue *res = graph.getResidueFromFile(comps[0], false);
            if (res != nullptr) {
                double score = atof(comps[1].c_str());
                scoreMap.setValue(res, score);
            }
        }
    }
    cout << "Finished loading scores" << endl;
    
    /*SeedGraph largestNeighborhood;
    Residue *residue;
    for (Residue *res: graph.getResidues()) {
        double *score = scoreMap.value(res);
        if (score != nullptr) {
            SeedGraph neighborhood = graph.neighborhood(res);
            if (neighborhood.residueSize() > 250) {
                largestNeighborhood = neighborhood;
                residue = res;
                break;
            }
        }
    }
    
    cout << "neighborhood to analyze has " << largestNeighborhood.residueSize() << " residues: " << graph.writeCodeForResidue(residue) << endl;*/
    
    PathFinder pathFinder(&graph, &scoreMap);
    cout << "Initialized" << endl;
    vector<pair<vector<Residue *>, mstreal>> paths = pathFinder.findPaths(20);
    cout << "Found " << paths.size() << " paths" << endl;
    
    ofstream ss(basePath + "bond_graphs/paths.txt", ios::out);
    if (!ss.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return 1;
    }
    
    for (auto pathInfo: paths) {
        vector<Residue *> path = pathInfo.first;
        //if (path.size() <= 3)
        //    continue;
        cout << pathInfo.second << ": ";
        ss << pathInfo.second << ",";
        for (Residue *res: path) {
            cout << *res << ",";
            ss << graph.writeCodeForResidue(res) << ",";
        }
        cout << endl;
        ss << endl;
    }
    
    

    return 0;
}
