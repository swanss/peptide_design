//
//  buildClusterTree.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 1/9/2020.
//

#include <stdio.h>
#include "clustertree.h"
#include "clusterutils.h"
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "fusehelper.h"
#include "seedgraph.h"
#include "seedscore.h"

using namespace std;

int main (int argc, char *argv[]) {
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Builds a cluster tree for fast fragment search retrieval.");
    opts.addOption("mode", "create = build a section of the whole tree, merge = combine smaller trees together", false);
    opts.addOption("pdbs", "A file containing absolute paths to all PDB structures to be used", false);
    opts.addOption("contacts", "A file containing absolute paths to contact files, corresponding exactly to the structures defined in the --pdbs parameter", false);

    opts.addOption("seeds", "A binary file containing seed structures. If provided, will produce a single-chain overlap tree with all fragments in a shared coordinate system.", false);
    opts.addOption("segLen", "Length of the segments to use for seed overlap clustering (default = 3).", false);
    opts.addOption("seedChain", "Name of the chain that represents the seed structure (default = '0').", false);

    opts.addOption("trees", "Directory in which trees to merge are stored. Required if mode = merge.", false);
    opts.addOption("out", "If no batching is specified, this is a file in which to write the output. Otherwise, the directory in which the output file should be stored. It will be named tree_{batch}.txt if mode = create, or tree_{mergeIndex}.txt if mode = merge.", false);
    opts.addOption("sizeLimit", "Number of fragments to add to the tree (only used if mode = create)", false);
    opts.addOption("batch", "The batch number (from 1 to numBatches), required if mode = create", false);
    opts.addOption("numBatches", "The number of batches (default 1)", false);
    opts.addOption("batchSplit", "Number of mini-batches that each worker should split the batch into (default 1).", false);
    opts.addOption("treeDim", "The maximum number of children of any node", true);
    opts.addOption("mergeSize", "The number of trees to merge", false);
    opts.addOption("mergeLimit", "The number of layers to cluster before performing a greedy reassignment of children", false);
    opts.addOption("mergeIndex", "The index of the merge job (required if mode = merge). Merges trees named tree_{(mergeIndex - 1) * treeDim}.txt through tree_{mergeIndex) * treeDim}.txt (i.e. it is 1-indexed).", false);
    opts.setOptions(argc, argv);
    
    string mode = opts.getString("mode", "create");
    string outPath = opts.getString("out");
    int treeDim = opts.getInt("treeDim");
    int batchSplit = opts.getInt("batchSplit", 1);
    
    FragmentFetcher *fetcher;
    StructuresBinaryFile *seedFile = nullptr;
    if (opts.isGiven("seeds")) {
        int segLen = opts.getInt("segLen", 3);
        string seedsPath = opts.getString("seeds");
        seedFile = new StructuresBinaryFile(seedsPath);
        seedFile->scanFilePositions();
        fetcher = new SingleFragmentFetcher(seedFile, segLen, opts.getString("seedChain", "0"), batchSplit);
    } else {
        MstUtils::assert(opts.isGiven("pdbs") && opts.isGiven("contacts"), "pdbs and contacts arguments must both be provided");
        string pdbsPath = opts.getString("pdbs");
        string contactsPath = opts.getString("contacts");
        fetcher = new PairFragmentFetcher(pdbsPath, contactsPath, 2, batchSplit);
    }

    if (mode == "create") {
        int batchIndex = opts.getInt("batch", 1) - 1;
        int numBatches = opts.getInt("numBatches", 1);
        if (batchIndex < 0 || numBatches < 1 || batchIndex >= numBatches) {
            cerr << "Batch index must be between 1 and numBatches" << endl;
            return 1;
        }
        cout << "Batch " << batchIndex << " of " << numBatches << endl;
        int sizeLimit = opts.getInt("sizeLimit", -1);

        for (int miniBatchIdx = 0; miniBatchIdx < batchSplit; miniBatchIdx++) {
            fetcher->beginNextSplit();
            if (batchSplit > 1)
                cout << "Mini-batch " << miniBatchIdx << " of " << numBatches * batchSplit << endl;
            BatchWorkerFragmentFetcher batchedFetcher(fetcher, numBatches * batchSplit, miniBatchIdx); 
            ClusterTree tree(&batchedFetcher, treeDim);
            tree.sizeLimit = sizeLimit;
            tree.cluster();
            if (opts.isGiven("numBatches") || opts.isGiven("batch") || opts.isGiven("batchSplit"))
                tree.write(MstSystemExtension::join(outPath, "tree_" + to_string(miniBatchIdx) + ".txt"));
            else
                tree.write(outPath);
        }

    } else if (mode == "merge") {
        string mergeTreesPath = opts.getString("trees", "");
        if (mergeTreesPath.size() == 0) {
            cerr << "Trees path must be provided when mode = merge" << endl;
            return 1;
        }
        int mergeIndex = opts.getInt("mergeIndex", 0) - 1;
        if (mergeIndex < 0) {
           cerr << "Merge index must be provided" << endl;
           return 1;
        } 
        int mergeSize = opts.getInt("mergeSize", treeDim);
        int mergeLimit = opts.getInt("mergeLimit", -1);
        
        vector<ClusterTree> trees;
        for (int i = mergeIndex * mergeSize; i < (mergeIndex + 1) * mergeSize; i++) {
            string treePath = MstSystemExtension::join(mergeTreesPath, "tree_" + to_string(i) + ".txt");
            if (!MstSys::fileExists(treePath)) {
                cerr << "File not found: " << treePath << endl;
                continue;
            }

            trees.emplace_back(fetcher, treeDim);
            trees.back().read(treePath);
            cout << "Read tree " << i << " of " << treeDim << endl;
        }
        
        if (trees.size() == 0) {
            cerr << "Nothing to merge" << endl;
            return 1;
        }
        ClusterTree mergedTree(trees, mergeLimit);
        if (opts.isGiven("numBatches") || opts.isGiven("batch") || opts.isGiven("batchSplit"))
            mergedTree.write(MstSystemExtension::join(outPath, "tree_" + to_string(mergeIndex) + ".txt"));
        else
            mergedTree.write(outPath);
    }
 
    /*// Make a query
    batchedFetcher->reset();
    if (!batchedFetcher->hasNext()) {
        cout << "No more results!" << endl;
        return 1;
    }
    auto item = batchedFetcher->next();
    ClusterSearchResults results = tree.search(item.first, 4.0);
    cout << "Found " << results.size() << " results" << endl;
    for (int i = 0; i < results.size(); i++) {
        cout << results.getFragmentInfo(i)->toString() << endl;
    }*/
    delete fetcher;
    if (seedFile != nullptr)
        delete seedFile;
    return 0;
}
