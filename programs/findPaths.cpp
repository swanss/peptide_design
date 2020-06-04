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
    opts.addOption("overlaps", "Path to a text file defining a cluster tree of overlaps", true);
    opts.addOption("out", "Path to a directory into which the fused seed path structures and scores will be written", true);
    opts.addOption("numPaths", "Number of paths to generate (linearly impacts running time - default 1000)", false);
    opts.addOption("ss", "Preferred secondary structure for paths (H, E, or O)", false);
    opts.setOptions(argc, argv);

    string targetPath = opts.getString("target");
    string binaryFilePath = opts.getString("seeds");
    string overlapTreePath = opts.getString("overlaps");
    string outputPath = opts.getString("out");

    if (!MstSys::fileExists(outputPath)) {
        MstSys::cmkdir(outputPath);
    }

    cout << "Loading seeds" << endl;
    StructuresBinaryFile seedFile(binaryFilePath);
    seedFile.scanFilePositions();
    SingleFragmentFetcher fetcher(&seedFile, 3, opts.getString("seedChain", "0"));
    ClusterTree overlapTree(&fetcher, 4, true); // 4 children per level; shared coordinate system
    overlapTree.read(overlapTreePath);

    Structure target(targetPath);

    // Remove native peptide from target
    if (opts.isGiven("peptideChain")) {
        Chain *peptide = target.getChainByID(opts.getString("peptideChain", "B"));
        if (peptide != nullptr)
            target.deleteChain(peptide);
    }

    // Scorer
    FragmentParams fParams(2, true);
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    StructureCompatibilityScorer scorer(&target, fParams, rParams, cParams, "/home/grigoryanlab/library/databases/dTERMen.databases/2019-01-22/dtermen.sim");

    int numPaths = opts.getInt("numPaths", 1000);

    ofstream out(MstSystemExtension::join(outputPath, "fused_paths.csv"), ios::out);
    if (!out.is_open())
        MstUtils::error("Could not open file stream");
    // CSV header
    out << "name,path,path_len,designability,num_contacts,num_designable,corroboration" << endl;

    // Stringent: 3-residue overlaps, 0.75A cutoff
    // Permissive: 3-residue overlaps, 1.25A cutoff
    PathSampler sampler(&target, &fetcher, &overlapTree, 3, 1.25, fetcher.getAllResidues());
    if (opts.isGiven("ss")) {
        sampler.preferredSecondaryStructure = new string(opts.getString("ss"));
    }

    int pathIndex = 0;
    StructuresBinaryFile fusedFile(MstSystemExtension::join(outputPath, "fused_paths.bin"), false);

    while (pathIndex < numPaths) {
        vector<PathResult> paths = sampler.sample(10);
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

            mstreal totalScore;
            int numContacts;
            int numDesignable;
            scorer.score(&fused, totalScore, numContacts, numDesignable);
            cout << "Score: " << totalScore << endl;
            out << totalScore << "," << numContacts << "," << numDesignable << "," << corroborationScore(path, 2) << endl;
            pathIndex++;
        }
    }

    cout << "Done" << endl;
    out.close();
    
    return 0;
}
