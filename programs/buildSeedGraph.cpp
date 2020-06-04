//
//  buildSeedGraph.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 1/28/19.
//

#include <stdio.h>
#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "clustertree.h"
#include "clusterutils.h"
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {

    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Builds connectivity graphs using a provided set of overlaps, and deposits them into sets of clusters ('subgraphs') and remaining 'chunks'.");
    opts.addOption("overlaps", "Directory defining overlaps (must be provided if --tree is not)", false);
    opts.addOption("tree", "File representing a cluster tree of k-mers (must be provided if --overlaps is not)", false);
    opts.addOption("bin", "Binary file containing seed structures", true);
    opts.addOption("out", "Directory into which to write seed graphs", true);
    opts.addOption("overlapSize", "The number of residues required to define an overlap (required if --tree is provided)", false);
    opts.addOption("overlapRMSD", "Maximum RMSD between sets of residues to be considered overlapping (default 1.0, only used if --tree is provided)", false);
    opts.addOption("minCosAngle", "Minimum cosine angle between Ca vectors allowed in a potential overlap (default -1.0, only used if --tree is provided)", false);
    opts.addOption("chunkSize", "Number of residues required to write out a chunk file (default is 100 if adj = same, 15 otherwise)", false);
    opts.addOption("adj", "If 'same' (default) then write graphs where adjacencies correspond to equivalent residues; if 'bond', write graphs where adjacencies are potential bonds", false);
    opts.setOptions(argc, argv);
    
    MstUtils::assert(opts.isGiven("overlaps") || opts.isGiven("tree"), "Either --overlaps or --tree must be provided");
    string overlapsPath = opts.getString("overlaps", "");
    string treePath = opts.getString("tree", "");
    string binaryPath = opts.getString("bin");
    string outPath = opts.getString("out");
    if (!MstSys::fileExists(outPath))
        MstSys::cmkdir(outPath);

    bool adjSameResidues = opts.getString("adj", "same") == "same";
    int chunkSize = opts.getInt("chunkSize", adjSameResidues ? 100 : 15);
    
    StructuresBinaryFile binaryFile(binaryPath);
    binaryFile.scanFilePositions();
    StructureCache *cache = new StructureCache(&binaryFile);
    cout << "Preloading structures..." << endl;
    cache->preloadFromBinaryFile();
    cout << "Done preloading structures" << endl;
    SeedGraph graph(false, cache);

    if (overlapsPath.size() > 0) {
        // Loading overlaps from explicit overlap paths
        int i = 0;
        string path = MstSystemExtension::join(overlapsPath, "overlaps" + to_string(i++) + ".txt");
        while (MstSystemExtension::fileExists(path)) {
            cout << "Loading candidates from batch " << i << endl;
            FuseCandidateFile file(path);
            graph.load(file); 
            path = MstSystemExtension::join(overlapsPath, "overlaps" + to_string(i++) + ".txt");
        }
    } else {
        // Loading overlaps from a cluster tree - any pair of nodes in the tree could be an overlap
        MstUtils::assert(opts.isGiven("overlapSize"), "--overlapSize must be specified if using a cluster tree");
        int overlapSize = opts.getInt("overlapSize");
        SingleFragmentFetcher fetcher(&binaryFile, cache, overlapSize, "0");
        ClusterTree tree(&fetcher, 4, true); 
        tree.read(treePath);
        cout << "Sorting subtrees by radius..." << endl;
        tree.sortSubtreesByRadius();
        cout << "Done sorting" << endl;

        graph.load(tree, overlapSize, opts.getReal("overlapRMSD", 1.0), opts.getReal("minCosAngle", -1.0));
    }
    cout << "Graph currently has " << graph.seedSize() << " seeds, now loading all seeds" << endl;
    graph.load(&binaryFile);
    
    cout << "Final graph has " << graph.seedSize() << " seeds and " << graph.residueSize() << " residues" << endl;
    SeedGraph graphToCluster = adjSameResidues ? graph.withAdjSameResidue() : graph;
    graphToCluster.write(MstSystemExtension::join(outPath, "whole_graph.txt"));
    
    cout << "Wrote whole graph, now finding subgraphs" << endl;
    vector<SeedGraph> subgraphs = graphToCluster.subgraphs();
    cout << "Found " << subgraphs.size() << " neighborhoods" << endl;

    // Test seed graph map
    /*unordered_set<Residue*> residues = subgraphs[0].getResidues();
    graph.setValue(*residues.begin(), "hello");
    for (Residue *res: residues) {
        string *val = graph.value(res);
        if (val != nullptr && *val != "hello")
            return 1;
    }
    cout << "All's well" << endl;*/
    
    // Write out neighborhoods of >= 15 residues; write out remainder in 15-cluster chunks
    SeedGraph currentChunk(false, graph.getStructures());
    int chunkIndex = 0;
    int chunkNumber = 0;
    
    int subgraphIndex = 0;
    vector<string> fileNames;
    
    for (SeedGraph g: subgraphs) {
        if (g.residueSize() < chunkSize) {
            currentChunk = currentChunk.unionWith(g);
            if (++chunkIndex == chunkSize) {
                string fileName = "chunk_" + to_string(chunkNumber++) + ".txt";
                currentChunk.write(MstSystemExtension::join(outPath, fileName));
                fileNames.push_back(fileName);
                currentChunk = SeedGraph(false, graph.getStructures());
                chunkIndex = 0;
            }
        } else {
            string fileName = "subgraph_" + to_string(subgraphIndex++) + ".txt";
            g.write(MstSystemExtension::join(outPath, fileName));
            fileNames.push_back(fileName);
        }
    }
    
    if (currentChunk.residueSize() > 0) {
        string fileName = "chunk_" + to_string(chunkNumber++) + ".txt";
        currentChunk.write(MstSystemExtension::join(outPath, fileName));
        fileNames.push_back(fileName);
    }
    
    // Write out tasks
    ofstream tasks(MstSystemExtension::join(outPath, "tasks.txt"), ios::out);
    if (!tasks.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return 1;
    }
    
    for (string fileName: fileNames) {
        tasks << fileName << endl;
    }
    
    cout << "done" << endl;
    delete cache;
    return 0;
}
