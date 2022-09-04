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
    opts.setTitle("Builds a graph describing seed residues and their potential connections and deposits them into sets of clusters ('subgraphs') and remaining 'chunks'.");
    opts.addOption("seedBin", "Binary file containing seed structures", true);
    opts.addOption("out", "File path to write out adjacency list for graph", true);
    opts.addOption("omitSeedsWithoutOverlaps","If provided, will omit seeds that do not overlap any other seed.",false);
    opts.addOption("overlaps", "The path to a directory (e.g. path/to/dir where the directory contains N files named with ascending values overlaps1.csv,overlaps2.csv,...,overlapsN.csv) or the path to a file defining overlaps (can have any name).", false);
    opts.addOption("tree", "File representing a cluster tree of k-mers (must be provided if --overlaps is not)", false);
    opts.addOption("adj", "If 'bond' (default) then write graphs where adjacencies are potential bonds; if 'same', write graphs where adjacencies correspond to equivalent residues", false);
    opts.addOption("subgraphDir", "If provided, write any disjoint subgraphs with at least --subgraphSize residues to this directory", false);
    opts.addOption("subgraphSize", "Number of residues required to write out a subgraph (default 100 for adj = same, 15 for adj = bond", false);
    opts.addOption("overlapSize", "The number of residues required to define an overlap (required if --tree is provided)", false);
    opts.addOption("overlapRMSD", "Maximum RMSD between sets of residues to be considered overlapping (default 1.0, only used if --tree is provided)", false);
    opts.addOption("minCosAngle", "Minimum cosine angle between Ca vectors allowed in a potential overlap (default -1.0, only used if --tree is provided)", false);
    opts.setOptions(argc, argv);
    
    MstUtils::assert(opts.isGiven("overlaps") || opts.isGiven("tree"), "Either --overlaps or --tree must be provided");
    string overlapsPath = opts.getString("overlaps", "");
    string treePath = opts.getString("tree", "");
    string binaryPath = opts.getString("seedBin");
    string outPath = opts.getString("out");
    bool omitSeedsWithoutOverlaps = opts.isGiven("omitSeedsWithoutOverlaps");
    string subgraphPath = opts.getString("subgraphDir", "");
    if (!subgraphPath.empty() && !MstSys::fileExists(subgraphPath))
        MstSys::cmkdir(subgraphPath);

    bool adjSameResidues = opts.getString("adj", "bond") == "same";
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
            cout << "Attempting to load overlaps from directory: " << overlapsPath << endl;
            int i = 1;
            string path = MstSystemExtension::join(overlapsPath, "overlaps" + to_string(i) + ".csv");
            while (MstSystemExtension::fileExists(path)) {
                FuseCandidateFile file(path);
                graph.load(file);
                i++;
                path = MstSystemExtension::join(overlapsPath, "overlaps" + to_string(i) + ".csv");
            }
            if (i==1) MstUtils::error("Expected to find at least one file containing overlaps with the name: "+path,"buildSeedGraph::main");
        } else {
            // All overlaps stored in one file
            cout << "Loading overlaps from file: " << overlapsPath << endl;
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
    if (graph.edgeSize() == 0) MstUtils::error("No seed overlaps found","buildSeedGraph::main"); 
    cout << "Graph currently has " << graph.seedSize() << " seeds, " << graph.residueSize() << " residues, and " << graph.edgeSize() << " edges." << endl;

    if (omitSeedsWithoutOverlaps) {
        cout << "Omitting seeds without overlaps to other seeds" << endl;
    } else {
        cout << "Now adding seeds that do not overlap other seeds" << endl;
        graph.loadCache();
    }
    
    cout << "Final graph has " << graph.seedSize() << " seeds, " << graph.residueSize() << " residues, and " << graph.edgeSize() << " edges" << endl;
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
    
    cout << "Done!" << endl;
    delete cache;
    return 0;
}
