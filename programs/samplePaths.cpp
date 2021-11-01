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
        if (coveredResidues.count(code) != 0) {
            continue;
        }
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
    opts.setTitle("Samples random paths from a seed graph/cluster tree and fuses the residues together into a peptide backbone structure.");
    opts.addOption("targetPDB", "Path to a PDB file of the target protein", true);
    opts.addOption("peptideChainID", "Single letter peptide chain ID. Will remove the peptide before generating seeds. Only necessary if the target_pdb contains a peptide chain that should be removed", false);
    opts.addOption("seedBin", "Path to a binary file containing seed structures", true);
    opts.addOption("seedChain", "Chain ID for the seed structures (default '0')", false);
    opts.addOption("seedGraph", "Path to a text file defining a seed graph", false);
    opts.addOption("numPaths", "Number of paths to generate. (default 200)", false);
    opts.addOption("minLength", "The minimum residue length of the sampled paths. (default 15)", false);
    opts.addOption("reqSeed", "The name of a seed in the binary file that all paths should extend from (optional)",false);
    opts.addOption("reqSeedSel", "A selection that specifies the residues in reqSeed that should always be included in sampled paths. Must be a continuous range: e.g. resid 3-5. (note: 'chain 0' is always assumed)",false);
    opts.addOption("fixedSeed", "If residues from the specified seed are included in a path, they will be fixed during fusing.",false);
    opts.addOption("numResFlank", "The number of flanking residues taken on each side of a contact when defining fragments during scoring (default 2)",false);
    opts.addOption("ss", "Preferred secondary structure for paths (H, E, or O)", false);
    opts.addOption("scorePaths", "Instead of sampling new paths from the graph, scores pre-defined paths. path format: seed_A:residue_i;seed_B:residue_j;etc...", false);
    opts.addOption("scoreStructures", "Instead of sampling new paths from the graph, loads structures, and scores.", false);
    opts.addOption("writeTopology","If provided, will write out the complete seed topology of each path. Useful for debugging.",false);
    opts.addOption("config", "The path to a configfile",true);
    opts.addOption("base", "Prepended to filenames",true);
    opts.setOptions(argc, argv);

    if (opts.isGiven("scoreStructures") && opts.isGiven("writeTopology")) MstUtils::error("Warning: will not write out topology in score structures mode");
    
    string targetPath = opts.getString("targetPDB");
    string binaryFilePath = opts.getString("seedBin");
    int numPaths = opts.getInt("numPaths", 200);
    int minLength = opts.getInt("minLength", 15);
    string configFilePath = opts.getString("config");
    string base = opts.getString("base");
    string seedChain = opts.getString("seedChain", "0");
    int flankingRes = opts.getInt("numResFlank",2);
    bool writeTopology = opts.isGiven("writeTopology");
    
    //The base name sets the seed, since this varies between batches, this should give unique sampling
    hash<string> hash;
    int name_hash = hash(base);
    cout << "base name: " << base << " with hash: " << name_hash << endl;
    srand(name_hash);
    
    string outputPath = "./path_structures";
    if (!MstSys::fileExists(outputPath)) {
        MstSys::cmkdir(outputPath);
    }
    string seedOutputPath = "./topology_seed_segments";
    if (writeTopology && !MstSys::fileExists(seedOutputPath)) {
        MstSys::cmkdir(seedOutputPath);
    }

    Structure target(targetPath);

    // Remove native peptide from target
    if (opts.isGiven("peptideChainID")) {
        Chain *peptide = target.getChainByID(opts.getString("peptideChainID", "B"));
        if (peptide != nullptr)
            target.deleteChain(peptide);
    }

    configFile config(configFilePath);
    RotamerLibrary RL(config.getRL());
    contactParams cParams(3.5,0.01,0.01);
    contactCounter cCounter(&target,&RL,cParams,0);
    
    StructuresBinaryFile* seedFile = nullptr;
    SeedGraphPathSampler* sampler = nullptr;
    StructureCache* cache = nullptr;
    SeedGraph* seedG = nullptr;
    
    if (opts.isGiven("seedGraph") && !opts.isGiven("scoreStructures")) {
        cout << "Loading seeds" << endl;
        seedFile = new StructuresBinaryFile(binaryFilePath);
        seedFile->scanFilePositions();
        
        cout << "Loading graph.." << endl;
        cache = new StructureCache(seedFile);
        seedG = new SeedGraph(opts.getString("seedGraph"), false, cache);
        
        int overlapLength = 1;
        sampler = new SeedGraphPathSampler(&target,seedG,overlapLength);
        if (opts.isGiven("ss")) {
            sampler->preferredSecondaryStructure = new string(opts.getString("ss"));
        }
        if (opts.isGiven("reqSeed")) {
            string reqSeedName = opts.getString("reqSeed");
            Structure* reqSeed = cache->getStructure(reqSeedName);
            vector<Residue*> residues;
            if (opts.isGiven("reqSeedSel")) {
                string reqSeedSel = opts.getString("reqSeedSel");
                selector sel(*reqSeed);
                residues = sel.selectRes(reqSeedSel);
                cout << "Select residues: ";
                for (Residue* R : residues) {
                    cout << R->getChainID() << R->getNum() << " ";
                }
                cout << endl;
            } else {
                Chain* C = reqSeed->getChainByID(seedChain);
                residues = C->getResidues();
            }
            sampler->setStartingPathResidues(residues);
        }
        if (opts.isGiven("fixedSeed")) {
            string fixedSeedName = opts.getString("fixedSeed");
            sampler->addFixedSeed(fixedSeedName);
        }
        sampler->setMinimumLength(minLength);
    }
    
    // Sample paths
    ofstream out(base+"_fused_paths.csv", ios::out);
    if (!out.is_open())
        MstUtils::error("Could not open file stream");
    // CSV header
    out << "name,path,path_len,fuser_score,rmsd_score,total_rmsd_score,bond_score,angle_score,dihedral_score,interchain_clash,intrachain_clash,num_interface_contacts,corroboration" << endl;
    
    cout << "Sample, fuse, and score " << numPaths << " fused paths..." << endl;
    int pathIndex = 0;
    while (pathIndex < numPaths) {
        vector<PathResult> paths;
        if (opts.isGiven("scorePaths")) {
            vector<string> path_strings = MstUtils::fileToArray(opts.getString("scorePaths"));
            paths = sampler->fusePaths(path_strings);
            numPaths = paths.size();
        } else if (opts.isGiven("scoreStructures")) {
            vector<string> structure_paths = MstUtils::fileToArray(opts.getString("scoreStructures"));
            vector<Structure> structures;
            for (string structure_path : structure_paths) {
                Structure S(structure_path);
                S.setName(MstSystemExtension::fileName(S.getName()));
                //verify that there is a single chain with ID = 0
                if (S.chainSize() != 1) MstUtils::error("Structures provided for scoring should have a single chain");
                if (S.getChainByID("0") == NULL) MstUtils::error("Structures provided for scoring must have a chain with ID = 0");
                vector<Residue*> empty; vector<vector<Residue*>> empty2;
                paths.emplace_back(empty,S,empty2,0,fusionOutput(),0,0);
            }
            numPaths = paths.size();
        } else {
            paths = sampler->sample(min(100,numPaths-pathIndex));
        }
        for (PathResult &path_result: paths) {
            cout << "Path: " << pathIndex << endl;
            string name = base + "_fused-path_" + to_string(pathIndex);
            string name_whole = base + "_fused-path-and-context_" + to_string(pathIndex);
            if (opts.isGiven("scoreStructures")) out << path_result.getFusedStructure().getName() << ",";
            else out << name << ",";
            
            for (Residue *res: path_result.getOriginalResidues()) {
                out << res->getStructure()->getName() << ":" << res->getResidueIndex() << ";";
            }
            out << "," << path_result.size() << ",";
            
            //report the fusion score (and its individual components)
            fusionOutput fuserScore = path_result.getFuserScore();
            out << fuserScore.getScore() << "," << fuserScore.getRMSDScore() << "," << fuserScore.getTotRMSDScore() << ",";
            out << fuserScore.getBondScore() << "," << fuserScore.getAngleScore() << "," << fuserScore.getDihedralScore() << ",";
            
            //report clashes
            out << path_result.getInterchainClash() << "," << path_result.getIntrachainClash() << ",";
            
            //get the path_and_context structure (including target-aligned residues) and path only
            Structure& path_and_context = path_result.getFusedStructure();
            Structure path_only;
            path_result.getFusedPathOnly(path_only);
            path_only.setName(name);
            path_and_context.setName(name_whole);

            int numContacts = 0;
            numContacts = cCounter.countContacts(&path_only);
            out << numContacts << ",";
            
            out << corroborationScore(path_result, 2) << endl;
            
            // Write out the PDBs
            path_and_context.writePDB(MstSystemExtension::join(outputPath,path_and_context.getName()+".pdb"));
            path_only.writePDB(MstSystemExtension::join(outputPath,path_only.getName()+".pdb"));
            
            if (writeTopology) {
                string seedOutputPathId = seedOutputPath + "/path_" + MstUtils::toString(pathIndex);
                if (!MstSys::fileExists(seedOutputPathId)) {
                    MstSys::cmkdir(seedOutputPathId);
                }
                vector<Structure> seedPathSegments = path_result.getTopologySeedResidues();
                for (Structure& S : seedPathSegments) {
                    S.writePDB(seedOutputPathId + "/" + S.getName()+".pdb");
                }
            }
            
            pathIndex++;
            if (pathIndex >= numPaths) break;
        }
    }
    
    if (!opts.isGiven("scoreStructures")) sampler->reportSamplingStatistics();
    
    if (opts.isGiven("seedGraph")) {
        delete cache;
        delete seedG;
        delete sampler;
    }

    cout << "Done" << endl;
    out.close();
    
    return 0;
}
