//
//  testClusterTree.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 1/6/2020.
//

#include <stdio.h>
#include <chrono>
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
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    //string binaryFilePath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seeds/1DUY_10000_200117/termextension_output/extendedfragments/extendedfragments.bin";    
    string pdbPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/cluster_db/nonredundantPDBs.list";
    string contactsPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/cluster_db/nonredundantContacts.list";
    int batchSize = atoi(argv[1]);
    int batchIndex = atoi(argv[2]) - 1;
    //string outPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/1DUY_10000_200117/tree.txt";

    /*StructuresBinaryFile binaryFile(binaryFilePath);
    binaryFile.scanFilePositions();
    SingleFragmentFetcher fetcher(&binaryFile, 3, "0");*/
    PairFragmentFetcher fetcher(pdbPath, contactsPath, 2);
    //FragmentFetcher *batchedFetcher = new BatchWorkerFragmentFetcher(&fetcher, batchSize, batchIndex); 
    
    // Test pairwise iteration
    ClusterTree tree(&fetcher, 4);
    tree.sizeLimit = 256;
    tree.cluster();
    cout << tree.toString() << endl;
    ClusterPairIterator it(tree);
    int visitCount = 0;
    while (it.hasNext()) {
        auto item = it.next();
        visitCount++;
    }
    cout << visitCount << " total visits" << endl;
    ClusterPairIterator it2(tree);
    visitCount = 0;
    while (it2.hasNext()) {
        auto item = it2.next();
        if (visitCount % 3 == 2)
            it2.skipSecondSubtree();
        visitCount++;
    }
    cout << visitCount << " total visits" << endl;

    // Merging
    /*vector<ClusterTree> trees;
    for (int i = 0; i < 4; i++) {
        trees.emplace_back(batchedFetcher, 4);
        trees.back().read("/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/cluster_db/test_tree_" + to_string(i) + ".txt");
        cout << "Read tree " << i << " of 4" << endl;
    }
    ClusterTree mergedTree(trees);
    mergedTree.write(outPath);*/

    /*ClusterTree tree(&fetcher, 4, true); // Shared coordinates = true
    tree.cluster();
    //tree.read(outPath);
    //cout << tree.toString() << endl;
    cout << "Writing" << endl;
    tree.write(outPath);*/

    // Make a query
    /*ClusterTree tree(batchedFetcher, 16);
    tree.read("/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/cluster_db/trees_ancestor_rmsd/trees_16/stage_3/tree_0.txt");
    batchedFetcher->reset();
    int numQueries = 0;    
    while (batchedFetcher->hasNext()) {
        if (rand() / (float)RAND_MAX > 0.1) {
            batchedFetcher->skip();
            continue;
        }
        auto item = batchedFetcher->next();

        for (double cutoff = 0.5; cutoff < 2.5; cutoff += 0.5) {
            cout << "Searching for fragment " << item.second->toString() << " with cutoff " << cutoff << endl;
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            ClusterSearchResults results = tree.search(item.first, cutoff);
            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            cout << duration_cast<seconds>( endTime - startTime ).count() << " seconds" << endl;
        }
        if (++numQueries >= 30)
            break;
    }*/
    // ClusterSearchResults results = tree.search(batchedFetcher->next().first, 1.5);
    // cout << "Found " << results.size() << " results" << endl;
    // if (!batchedFetcher->hasNext()) {
    //     cout << "No more results!" << endl;
    //     return 1;
    // }
    // ClusterSearchResults newResults = tree.search(batchedFetcher->next().first, 1.5);
    // cout << "Found " << newResults.size() << " results" << endl;
    // delete batchedFetcher;
    return 0;
}
