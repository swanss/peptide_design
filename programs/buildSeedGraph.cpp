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
    opts.addOption("overlaps", "Directory or file defining overlaps (must be provided if --tree is not)", false);
    opts.addOption("tree", "File representing a cluster tree of k-mers (must be provided if --overlaps is not)", false);
    opts.addOption("bin", "Binary file containing seed structures", true);
    opts.addOption("out", "File path to write out adjacency list for graph", true);
    opts.addOption("adj", "If 'same' (default) then write graphs where adjacencies correspond to equivalent residues; if 'bond', write graphs where adjacencies are potential bonds", false);
    opts.addOption("subgraphDir", "If provided, write any disjoint subgraphs with at least --subgraphSize residues to this directory", false);
    opts.addOption("subgraphSize", "Number of residues required to write out a subgraph (default 100 for adj = same, 15 for adj = bond", false);
    opts.addOption("overlapSize", "The number of residues required to define an overlap (required if --tree is provided)", false);
    opts.addOption("overlapRMSD", "Maximum RMSD between sets of residues to be considered overlapping (default 1.0, only used if --tree is provided)", false);
    opts.addOption("minCosAngle", "Minimum cosine angle between Ca vectors allowed in a potential overlap (default -1.0, only used if --tree is provided)", false);
    opts.setOptions(argc, argv);
    
    MstUtils::assert(opts.isGiven("overlaps") || opts.isGiven("tree"), "Either --overlaps or --tree must be provided");
    string overlapsPath = opts.getString("overlaps", "");
    string treePath = opts.getString("tree", "");
    string binaryPath = opts.getString("bin");
    string outPath = opts.getString("out");
    string subgraphPath = opts.getString("subgraphDir", "");
    if (!subgraphPath.empty() && !MstSys::fileExists(subgraphPath))
        MstSys::cmkdir(subgraphPath);

    bool adjSameResidues = opts.getString("adj", "same") == "same";
    int subgraphSize = opts.getInt("subgraphSize", adjSameResidues ? 100 : 15);
    
    StructuresBinaryFile binaryFile(binaryPath);
    binaryFile.scanFilePositions();
    StructureCache *cache = new StructureCache(&binaryFile);
    cout << "Preloading structures..." << endl;
    cache->preloadFromBinaryFile();
    cout << "Done preloading structures" << endl;
    SeedGraph graph(false, cache, true);

    if (overlapsPath.size() > 0) {
        // Loading overlaps from explicit overlap paths
        if (MstSys::isDir(overlapsPath)) {
            int i = 0;
            string path = MstSystemExtension::join(overlapsPath, "overlaps" + to_string(i++) + ".txt");
            while (MstSystemExtension::fileExists(path)) {
                cout << "Loading overlaps from batch " << i << endl;
                FuseCandidateFile file(path);
                graph.load(file);
                path = MstSystemExtension::join(overlapsPath, "overlaps" + to_string(i++) + ".txt");
            }
        } else {
            // All overlaps stored in one file
            cout << "Loading overlaps..." << endl;
            FuseCandidateFile file(overlapsPath);
            graph.load(file);
            cout << "Done loading overlaps" << endl;
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
    graphToCluster.write(outPath);
    
    if (!subgraphPath.empty()) {
        cout << "Wrote whole graph, now finding subgraphs" << endl;
        vector<SeedGraph> subgraphs = graphToCluster.subgraphs();
        cout << "Found " << subgraphs.size() << " subgraphs" << endl;
        
        // Write out neighborhoods of >= subgraphSize residues
        SeedGraph currentChunk(false, graph.getStructures());
        int chunkIndex = 0;
        int chunkNumber = 0;
        
        int subgraphIndex = 0;
        vector<string> fileNames;
        
        for (SeedGraph g: subgraphs) {
            if (g.residueSize() >= subgraphSize) {
                string fileName = "subgraph_" + to_string(subgraphIndex++) + ".txt";
                g.write(MstSystemExtension::join(subgraphPath, fileName));
                fileNames.push_back(fileName);
            }
        }
        
        cout << subgraphIndex << " of the subgraphs meet the size threshold" << endl;
        
        // Write out a list of the subgraphs
        ofstream subgraphList(MstSystemExtension::join(subgraphPath, "subgraphs.list"), ios::out);
        if (!subgraphList.is_open()) {
            cerr << "couldn't open out stream" << endl;
            return 1;
        }
        
        for (string fileName: fileNames) {
            subgraphList << fileName << endl;
        }
    }
    
    cout << "done" << endl;
    delete cache;
    return 0;
}
