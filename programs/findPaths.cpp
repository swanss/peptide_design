//
//  findPaths.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 4/22/19.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "structure_iter.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "pathsampler.h"
#include "findpaths.h"
#include "utilities.h"
#include <unordered_set>

using namespace std;
using namespace std::chrono;

/**
 * Helper function that computes the corroboration score of a given path. The
 * corroboration score is the number of unique target residues that were used
 * as sources for the seeds that gave rise to the residues in the path.
 *
 * @param path the path to score
 * @param nearbyThreshold the number of residues to either side of each target
 *  residue to ignore in the final count. For example, if nearbyThreshold is 2
 *  and the target residue A53 is discovered, residues A51-A55 will be prevented
 *  from counting toward the corroboration score.
 */
int corroborationScore(PathResult &path, int nearbyThreshold) {
    int score = 0;
    vector<pair<string, int>> targetResidues;
    for (Residue *res: path.getOriginalResidues()) {
        string seedName = res->getStructure()->getName();
        auto code = getTargetResidueCode(seedName);
        targetResidues.push_back(code);
    }

    // Put target residues in sorted order before eliminating close neighbors
    sort(targetResidues.begin(), targetResidues.end(), [](const pair<string, int> &p1, const pair<string, int> &p2) {
        if (p1.first.compare(p2.first) < 0)    
            return true;
        return p1.second < p2.second;
    });
    unordered_set<pair<string, int>, pair_hash> coveredResidues;
    for (auto code: targetResidues) {
        if (coveredResidues.count(code) != 0)
            continue;
        score++;
        for (int i = code.second - nearbyThreshold; i < code.second + nearbyThreshold + 1; i++) {
            coveredResidues.insert(make_pair(code.first, i));
        }
    }

    return score;
}

int main (int argc, char *argv[]) {
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Finds overlaps between all pairs of the given seed residues.");
    opts.addOption("target", "Path to the target PDB structure file", true);
    opts.addOption("peptideChain", "Chain ID for the peptide chain in the target, if one exists - it will be removed", false);
    opts.addOption("seeds", "Path to a binary file containing seed structures", true);
    opts.addOption("seedChain", "Chain ID for the seed structures (default is '0')", false);
    opts.addOption("overlaps", "Path to a text file defining a cluster tree of overlaps", false);
    opts.addOption("seedGraph", "Path to a text file defining a seed graph", false);
    opts.addOption("out", "Path to a directory into which the fused seed path structures and scores will be written", true);
    opts.addOption("numPaths", "Number of paths to generate (linearly impacts running time - default 1000)", false);
    opts.addOption("req_seed", "The name of a seed in the binary file that all paths should extend",false);
    opts.addOption("ss", "Preferred secondary structure for paths (H, E, or O)", false);
    opts.addOption("config", "The path to a configfile",true);
    opts.addOption("noScore", "If provided, disable designability and contact scoring", false);
    opts.setOptions(argc, argv);

    if (opts.isGiven("overlaps") == opts.isGiven("seedGraph")) MstUtils::error("Either 'overlaps' or 'seedGraph' must be provided, but not both.");
    
    string targetPath = opts.getString("target");
    string binaryFilePath = opts.getString("seeds");
    string overlapTreePath = opts.getString("overlaps");
    string outputPath = opts.getString("out");
    string configFilePath = opts.getString("config");
    
    string seedChain = opts.getString("seedChain", "0");
    bool shouldScore = !opts.isGiven("noScore");
    
    if (!MstSys::fileExists(outputPath)) {
        MstSys::cmkdir(outputPath);
    }

    Structure target(targetPath);

    // Remove native peptide from target
    if (opts.isGiven("peptideChain")) {
        Chain *peptide = target.getChainByID(opts.getString("peptideChain", "B"));
        if (peptide != nullptr)
            target.deleteChain(peptide);
    }

    // Scorer
    StructureCompatibilityScorer *scorer = nullptr;
    if (shouldScore) {
        FragmentParams fParams(2, true);
        rmsdParams rParams(1.2, 15, 1);
        contactParams cParams;
        scorer = new StructureCompatibilityScorer(&target, fParams, rParams, cParams, configFilePath);
    }

    int numPaths = opts.getInt("numPaths", 1000);
    
    cout << "Loading seeds" << endl;
    StructuresBinaryFile seedFile(binaryFilePath);
    seedFile.scanFilePositions();
    
    PathSampler* sampler;
    SingleFragmentFetcher* fetcher;
    ClusterTree* overlapTree;
    StructureCache* cache;
    SeedGraph* seedG;
    
    // Handle configuration options separately for each input type, since different
    // sets of options may be available
    if (opts.isGiven("overlaps")) {
        cout << "Loading cluster tree..." << endl;
        fetcher = new SingleFragmentFetcher(&seedFile, 3, seedChain);
        overlapTree = new ClusterTree(fetcher, 4, true); // 4 children per level; shared coordinate system
        overlapTree->read(overlapTreePath);
        
        // Stringent: 3-residue overlaps, 0.75A cutoff
        // Permissive: 3-residue overlaps, 1.25A cutoff
        ClusterTreePathSampler *cSampler = new ClusterTreePathSampler(&target, fetcher, overlapTree, 3, 1.25, fetcher->getAllResidues());
        if (opts.isGiven("ss")) {
            cSampler->preferredSecondaryStructure = new string(opts.getString("ss"));
        }
        if (opts.isGiven("req_seed")) {
            MstUtils::error("req_seed parameter not implemented for cluster tree path sampling");
        }
        
        sampler = cSampler;
        
    } else if (opts.isGiven("seedGraph")) {
        cout << "Loading graph.." << endl;
        cache = new StructureCache(&seedFile);
        seedG = new SeedGraph(false,cache);
        
        SeedGraphPathSampler *gSampler = new SeedGraphPathSampler(&target,seedG);
        if (opts.isGiven("ss")) {
            gSampler->preferredSecondaryStructure = new string(opts.getString("ss"));
        }
        if (opts.isGiven("req_seed")) {
            string reqSeedName = opts.getString("req_seed");
            Structure* reqSeed = cache->getStructure(reqSeedName);
            gSampler->setStartingSeed(reqSeed,seedChain);
        }

        sampler = gSampler;
    }
    
    ofstream out(MstSystemExtension::join(outputPath, "fused_paths.csv"), ios::out);
    if (!out.is_open())
        MstUtils::error("Could not open file stream");
    // CSV header
    out << "name,path,path_len,designability,num_contacts,num_designable,corroboration" << endl;
    
    int pathIndex = 0;
    StructuresBinaryFile fusedFile(MstSystemExtension::join(outputPath, "fused_paths.bin"), false);

    while (pathIndex < numPaths) {
        vector<PathResult> paths = sampler->sample(10);
        for (PathResult &path: paths) {
            string name = "fused_path_" + to_string(pathIndex);
            out << name << ",";
            for (Residue *res: path.getOriginalResidues()) {
                out << res->getStructure()->getName() << ":" << res->getResidueIndex() << ";";
            }
            out << "," << path.size() << ",";

            // Score the path
            Structure fused;
            path.getFusedPathOnly(fused);
            fused.setName(name);
            fusedFile.appendStructure(&fused);
            Structure &fullFused = path.getFusedStructure();

            // Write zeros for scores if scoring is disabled (speeds up performance)
            mstreal totalScore = 0.0;
            int numContacts = 0;
            int numDesignable = 0;
            if (shouldScore) {
                scorer->score(&fused, totalScore, numContacts, numDesignable);
                cout << "Score: " << totalScore << endl;
            }
            out << totalScore << "," << numContacts << "," << numDesignable << "," << corroborationScore(path, 2) << endl;
            pathIndex++;
        }
    }
    
    if (opts.isGiven("overlaps")) {
        delete fetcher;
        delete overlapTree;
    } else if (opts.isGiven("seedGraph")) {
        delete cache;
        delete seedG;
    }
    delete sampler;
    if (shouldScore)
        delete scorer;

    cout << "Done" << endl;
    out.close();
    
    return 0;
}
