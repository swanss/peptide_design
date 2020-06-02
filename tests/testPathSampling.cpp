//
//  testPathSampling.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 1/23/2020.
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
#include "pathsampler.h"
#include "Util.h"

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    
    string structureID = argv[1];
    string peptideChain = argv[2];
    string binaryFilePath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seeds/" + structureID + "_200204/termextension_output/extendedfragments.bin";    
    string pdbPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/cluster_db/nonredundantPDBs.list";
    string contactsPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/cluster_db/nonredundantContacts.list";
    string treePath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/cluster_db/trees_ancestor_rmsd/trees_16/stage_3/tree_0.txt";
    string overlapTreePath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/" + structureID + "_200204/tree.txt";
    string graphPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/" + structureID + "_200204/bond_graphs/whole_graph.txt";
    string outputPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/" + structureID + "_200204";
    string targetPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/seed_graphs/" + structureID + "_200204/" + structureID + ".pdb";

    /*PairFragmentFetcher fetcher(pdbPath, contactsPath, 2);
    FragmentFetcher *batchedFetcher = new BatchWorkerFragmentFetcher(&fetcher, batchSize, batchIndex); 

    ClusterTree tree(batchedFetcher, 16);
    tree.read(treePath);*/

    cout << "Loading seeds" << endl;
    StructuresBinaryFile seedFile(binaryFilePath);
    seedFile.scanFilePositions();
    SingleFragmentFetcher fetcher(&seedFile, 3, "0");
    ClusterTree overlapTree(&fetcher, 4, true);
    overlapTree.read(overlapTreePath);

    //StructureCache cache(&seedFile);
    //cout << "Loading graph" << endl;
    //SeedGraph graph(graphPath, false, &cache);
    Structure target(targetPath);

    // Remove native peptide from target
    Chain *peptide = target.getChainByID(peptideChain);
    if (peptide != nullptr)
        target.deleteChain(peptide);

    // Scorer
    FragmentParams fParams(2, true);
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    StructureCompatibilityScorer scorer(&target, fParams, rParams, cParams, "/home/grigoryanlab/library/databases/dTERMen.databases/2019-01-22/dtermen.sim");

    // Check for clashes with seeds
    /*seedFile.reset();
    int numClash = 0;
    int numTested = 0;
    while (seedFile.hasNext()) {
        Structure *s = seedFile.next();
        Structure seedOnly;
        seedOnly.appendChain(new Chain(*(s->getChainByID("0"))));
        if (scorer.clashes(&seedOnly)) {
            numClash++;
            cout << s->getName() << " clashes with target" << endl;
        }
        numTested++;
        delete s;
    }
    cout << numClash << " out of " << numTested << " clash with target" << endl;
    return 0;*/

    int numPaths = 250;

    ofstream out(MstSystemExtension::join(outputPath, "fused_paths_unique_seeds.list"), ios::out);
    if (!out.is_open())
        MstUtils::error("Could not open file stream");

    PathSampler sampler(&target, &fetcher, &overlapTree, 3, 0.75, fetcher.getAllResidues()); // &graph, &cache);
    sampler.uniqueSeeds = true;

    int pathIndex = 0;
    StructuresBinaryFile fusedFile(MstSystemExtension::join(outputPath, "fused_paths_unique_seeds.bin"), false);

    while (pathIndex < numPaths) {
        vector<PathResult> paths = sampler.sample(10);
        for (PathResult &path: paths) {
            string name = "fused_path_" + to_string(pathIndex);
            out << name << ",";
            for (Residue *res: path.getOriginalResidues()) {
                out << res->getStructure()->getName() << ":" << res->getResidueIndex() << ";";
            }
            out << ",";

            // Score the path
            Structure fused;
            path.getFusedPathOnly(fused);
            fused.setName(name);
            fusedFile.appendStructure(&fused);
            Structure &fullFused = path.getFusedStructure();

            mstreal totalScore;
            int numContacts;
            int numDesignable;
            scorer.score(&fused, totalScore, numContacts, numDesignable);
            cout << "Score: " << totalScore << endl;
            out << totalScore << "," << numContacts << "," << numDesignable << endl;
            pathIndex++;
        }
    }

    cout << "Done" << endl;
    out.close();
    return 0;
}
