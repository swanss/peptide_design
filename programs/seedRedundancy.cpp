//
//  seedRedundancy.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 2/21/19.
//

#include <stdio.h>
#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "fusehelper.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "findpaths.h"
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    string basePath = "/home/ifs-users/venkats/dtermen/binding_site_new_seeds/";
    string dataPath = "/home/ifsdata/scratch/grigoryanlab/swans/peptide_design_benchmark/seed_retrieval/target_centered_approach/retrieveSeedsTest3/1DUY/output/seeds/";
    string filesPath = basePath + "seeds_gen0.csv";
    
    // mode = 1: corresponding-residue graph; mode = 2: adjacent-residue graph
    int mode = 0;
    std::istringstream ss(argv[1]);
    if (!(ss >> mode)) {
        std::cerr << "Invalid number: " << argv[1] << '\n';
    } else if (!ss.eof()) {
        std::cerr << "Trailing characters after number: " << argv[1] << '\n';
    }
    
    string outPath(argv[2]);

    bool useSameResidues = (mode == 1);
    SeedGraph graph(useSameResidues);
    for (int i = 0; i < 20; i++) {
        cout << "Loading candidates from batch " << i << endl;
        FuseCandidateFile file(basePath + "overlaps_res2/results" + to_string(i) + ".txt");
        graph.load(file, dataPath);
    }
    
    cout << "Built graph with " << graph.seedSize() << " seeds and " << graph.residueSize() << " residues" << endl;
    unordered_set<Residue *> seenResidues;
    int numFlank = 2;
    ofstream outFile(outPath, ios::out);
    if (!outFile.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return 1;
    }

    for (Residue *res: graph.getResidues()) {
        if (seenResidues.count(res) != 0) continue;
        
        vector<Residue *> analogs = graph.similarContactResidues(res, numFlank, 1.0f);
        seenResidues.insert(analogs.begin(), analogs.end());
        seenResidues.insert(res);
        outFile << MstSystemExtension::fileName(res->getStructure()->getName()) << ":" << res->getName() << ":" << res->getResidueIndex() << ",";
        for (Residue *res2: analogs) {
            outFile << MstSystemExtension::fileName(res2->getStructure()->getName()) << ":" << res2->getName() << ":" << res2->getResidueIndex();
            outFile << ",";
        }
        outFile << endl;
    }
    //graph.write(outPath);
    
    /*cout << "Finding neighborhoods" << endl;
    vector<SeedGraph> subgraphs = graph.subgraphs();
    int subgraphIndex = 0;
    
    ofstream ss(outPath, ios::out);
    if (!ss.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return 1;
    }
    
    for (SeedGraph g: subgraphs) {
        if (g.residueSize() > 1)
            cout << "Cluster with " << g.residueSize() << " residues" << endl;
        for (Residue *res: g.getResidues()) {
            ss << MstSystemExtension::pathBase(res->getStructure()->getName()) << ":" << res->getName() << ":" << res->getResidueIndex();
            ss << ",";
        }
        ss << endl;
    }*/
    cout << "Done" << endl;
    
    return 0;
}
