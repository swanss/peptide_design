//
//  sandbox.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 11/13/18.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstcondeg.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "fusehelper.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "ClusterSearch.h"
#include "Fragments.h"
#include <fstream>
#include <sstream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main (int argc, char *argv[]) {
    cout << "hello world!" << endl;
//    
//    // Read scoring mode from argv
//    int mode = 0;
//    istringstream ss(argv[1]);
//    if (!(ss >> mode)) {
//        std::cerr << "Invalid number: " << argv[1] << '\n';
//    } else if (!ss.eof()) {
//        std::cerr << "Trailing characters after number: " << argv[1] << '\n';
//    }
//    int numBatches = 0;
//    ss = istringstream(argv[2]);
//    if (!(ss >> numBatches)) {
//        std::cerr << "Invalid number: " << argv[2] << '\n';
//    } else if (!ss.eof()) {
//        std::cerr << "Trailing characters after number: " << argv[2] << '\n';
//    }
//    
//    int numResOverlap = 2;
//    string basePath = "/home/ifs-users/venkats/dtermen/binding_site_new_seeds/";
//    string dataPath = "/home/ifsdata/scratch/grigoryanlab/swans/peptide_design_benchmark/seed_retrieval/target_centered_approach/retrieveSeedsTest3/1DUY/output/seeds/"; //basePath + "1DUY_decorate_sample/";
//    string outPath = basePath + "seed_scores_no_clash/"; //"/home/ifs-users/venkats/dtermen/nearby_seeds_test/new_results_res" + to_string(numResOverlap) + "/";
//    string filesPath = basePath + "seeds_gen0.csv"; //"/home/ifsdata/scratch/grigoryanlab/swans/peptide_design_benchmark/seed_retrieval/library_centered_approach/1DUY_100k_terms-080318/1DUY_decorate_output/";
//    string fasstDB = "/home/grigoryanlab/library/databases/dTERMen.databases/2019-01-22/dtermen.sim";
//
//    //string seedPath = basePath + "nchains1_nres10_score-0.202_-0.200_rmsd0.198_term049253_chidsA_match23.pdb";
//    Structure target("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb");
//    removeSideChains(target);
//    cout << "Reading database" << endl;
//    ClusterDatabase cd("/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/cluster_db/clusterDBNoCentroid", 2); 
//    cout << "Read database" << endl;
//   
//    Structure testS = cd.getStructure(1, 1);
//    cout << "Test structure has " << testS.residueSize() << " residues" << endl;
// 
//    ConFind confind("/home/grigoryanlab/library/MST/testfiles/rotlib.bin", target);
//    contactList cl = confind.getContacts(target, 0.03, nullptr);
//    FragmentParams params(2, 2, true, true, true, false);
//    Fragmenter fragmenter(target, params);
//    fragmenter.fragment(cl);
//    vector<Structure> fragments = fragmenter.getFragmentStructures();
//    for (Structure s: fragments) {
//        cout << s.atomSize() << ", " << cd.numAtoms() << endl;
//        if (s.atomSize() == cd.numAtoms()) {
//            ClusterSearch search(&cd);
//            search.setQuery(s);
//            search.search(1.5);
//            cout << "Found " << search.numMatches() << " matches" << endl;
//            break;
//        }
//    }

    // 1/28/19: Score the results of fusing a subset of samples
    /*FuseCandidateFile file("/home/ifs-users/venkats/dtermen/nearby_seeds_test/subgraphs_binding_site_strict/used_candidates_small_sample.csv");
    FuseHelper helper;
    Structure *target = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb");
    FragmentParams fParams(2, false);
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    string fasstDB = "/home/grigoryanlab/library/databases/FASST/db-12657.bin";
    SeedScorer *scorer = new SequenceCompatibilityScorer(target, fParams, rParams, cParams, fasstDB, 0.4, 0.05, 0.25, 1, 8000, 0.7); // StructureCompatibilityScorer(target, fParams, rParams, cParams, fasstDB);

    outPath = "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/fuse_binding_site_seq_score/";
    SeedListFile summaryFile(outPath + "fuse_results.csv");
    summaryFile.setMetadataFields({ "pre_seq_score_1", "pre_seq_score_2", "post_seq_score", "total_score", "bond_score", "angle_score", "dihedral_score", "rmsd_score", "total_rmsd_score", "source_1", "source_chain_1", "source_2", "source_chain_2" });
    
    vector<FuseCandidate> candidateBatch;
    int fuseIndex = 0;
    while ((candidateBatch = file.read(1)).size() > 0) {
        FuseCandidate candidate = candidateBatch[0];
        candidate.loadStructures(dataPath);
        Structure scoringStruct1;
        Structure scoringStruct2;
        extractChains(*candidate.structure1, candidate.chain1, scoringStruct1);
        extractChains(*candidate.structure2, candidate.chain2, scoringStruct2);
        auto preScore1 = scorer->score(&scoringStruct1);
        auto preScore2 = scorer->score(&scoringStruct2);

        string pdbName = "fuse_result_" + to_string(fuseIndex++) + ".pdb";
        // Fused structure, chain IDs of seed, and fusion output indicating scores
        tuple<Structure, string, fusionOutput> result = helper.performFuse(candidate);
        get<0>(result).writePDB(outPath + pdbName);
        fusionOutput output = get<2>(result);
        Structure scoringStructFull;
        extractChains(get<0>(result), get<1>(result), scoringStructFull);
        auto postScore = scorer->score(&scoringStructFull);
        
        summaryFile.write(pdbName, get<1>(result), {
            // Scores to write along with the seed info
            writePerResidueScores(preScore1),
            writePerResidueScores(preScore2),
            writePerResidueScores(postScore),
            to_string(output.getScore()),
            to_string(output.getBondScore()),
            to_string(output.getAngleScore()),
            to_string(output.getDihedralScore()),
            to_string(output.getRMSDScore()),
            to_string(output.getTotRMSDScore()),
            candidate.file1,
            candidate.chain1,
            candidate.file2,
            candidate.chain2
        });
        candidate.freeStructures();
    }
    cout << "Fused " << fuseIndex << " candidates" << endl;
    delete target;
    delete scorer;*/
    
    // 3/18/19: Get contact scores
    /*Structure *target = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb");
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceCompatibilityScorer scorer(target, rParams, cParams, fasstDB, 2, 0, 0.4, 0.05, 0.25, 1, 8000, 0.7);
    scorer.noQueries = true;
    
    string graphsPath = basePath + "graphs/";
    vector<string> tasks = MstUtils::fileToArray(graphsPath + "tasks.txt");
    
    int i = 0;
    unordered_set<string> scoredSeeds;
    for (int taskIdx = mode - 1; taskIdx < tasks.size(); taskIdx += numBatches) {
        string graphFile = tasks[taskIdx];
        cout << "graph file: " << graphFile << endl;
        // Load adjacency graph
        SeedGraphMap<mstreal> graph;
        graph.read(MstSystemExtension::join(graphsPath, graphFile), dataPath);
        scorer.adjacencyGraph = &graph;
        StructureCache *structures = graph.getStructures();
        
        for (auto it = structures->begin(); it != structures->end(); it++) {
            Structure *seed = *it;
            if (scoredSeeds.count(seed->getName()) != 0) {
                cout << "redundant seed!" << endl;
                continue;
            }
            scoredSeeds.insert(seed->getName());
            
            if ((i++) % 100 == 0)
                cout << "Structure " << i << endl;
            if (seed->chainSize() > 1)
                continue;
            cout << "Loaded structure " << seed->getName() << endl;
            
            // Collect contacts from the seed
            scorer.collectContacts(seed);
        }
    }
    
    cout << "Done counting contacts!" << endl;
    scorer.writeContactCounts(outPath + "contact_scores/contact_scores_" + to_string(mode) + ".csv");
    
    delete target;*/

    // 3/9/19: Score all seeds in parallel
    /*Structure *target = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb");
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceCompatibilityScorer scorer(target, rParams, cParams, fasstDB, 2, 0, 0.4, 0.05, 0.25, 1, 8000, 0.7);
    string scoreWritePath = outPath + "all_fragment_scores/frag_scores_" + to_string(mode) + ".csv";
    scorer.scoresWritePath = &scoreWritePath;
    
    // Read contact counts
    string contactCountsPath = outPath + "contact_scores/"; //basePath + "seed_scores/contact_scores/";
    for (int i = 1; i <= 35; i++) {
        scorer.readContactCounts(contactCountsPath + "contact_scores_" + to_string(i) + ".csv");
    }

    string graphsPath = basePath + "graphs/";
    vector<string> tasks = MstUtils::fileToArray(graphsPath + "tasks.txt");
    
    double totalTime = 0.0;
    int numTimes = 0;
    
    ofstream outputSS(outPath + "seeds_scored_flank0_" + to_string(mode) + ".csv", ios::out);
    if (!outputSS.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return 1;
    }

    int i = 0;
    unordered_set<string> scoredSeeds;
    for (int taskIdx = mode - 1; taskIdx < tasks.size(); taskIdx += numBatches) {
        string graphFile = tasks[taskIdx];
        cout << "graph file: " << graphFile << endl;
        // Load adjacency graph
        SeedGraphMap<mstreal> graph;
        graph.read(MstSystemExtension::join(graphsPath, graphFile), dataPath);
        scorer.adjacencyGraph = &graph;
        StructureCache *structures = graph.getStructures();

        for (auto it = structures->begin(); it != structures->end(); it++) {
            Structure *seed = *it;
            if (scoredSeeds.count(seed->getName()) != 0) {
                cout << "redundant seed!" << endl;
                continue;
            }
            scoredSeeds.insert(seed->getName());
            
            if ((i++) % 100 == 0)
                cout << "Structure " << i << endl;
            if (seed->chainSize() > 1)
                continue;
            cout << "Loaded structure " << seed->getName() << endl;
            
            // Score the seed
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            auto result = scorer.score(seed);
            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            double time = duration_cast<seconds>( endTime - startTime ).count();
            cout << "time: " << time << endl;
            totalTime += time;
            numTimes += 1;
            
            //outputSeedFile.write(seed->getName(), seed->getChain(0).getID(), { writePerResidueScores(result) });
            for (auto resScore: result) {
                outputSS << graph.writeCodeForResidue(resScore.first) << "," << resScore.second << endl;
            }
        }
    }
    
    cout << "Done! Average time per seed: " << totalTime << "/" << numTimes << " = " << totalTime / numTimes << endl;
    cout << "Scoring statistics: " << endl;
    cout << "\t" << scorer.numResiduesScored() << " residues scored" << endl;
    cout << "\t" << scorer.numUniqueResiduesScored() << " unique residues scored" << endl;
    cout << "\t" << scorer.numSeedQueries() << " seed queries" << endl;
    cout << "\t" << scorer.numTargetQueries() << " target queries" << endl;
    
    delete target;*/
    
    // 1/25/19: Score all seeds
    /*Structure *target = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb");
    Structure *realPeptide = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/test_pose_4U6Y.pdb");
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceCompatibilityScorer scorer(target, rParams, cParams, fasstDB, 2, 0, 0.4, 0.05, 0.25, 1, 8000, 0.7);
    
    // Load adjacency graph
    SeedGraphMap<mstreal> graph;
    for (int i = 0; i < 20; i++) {
        cout << "Loading candidates from batch " << i << endl;
        FuseCandidateFile file(basePath + "overlaps_res2/results" + to_string(i) + ".txt");
        graph.load(file, dataPath);
    }
    scorer.adjacencyGraph = &graph;
    
    cout << "Loaded scorer" << endl;
    double totalTime = 0.0;
    int numTimes = 0;
    
    SeedListFile inputSeedFile(filesPath);
    SeedListFile outputSeedFile(basePath + "seeds_scored_flank0_" + to_string(mode) + ".csv");
    outputSeedFile.setMetadataFields({ mode == 1 ? "sequence_score" : "structure_score" });
    auto seedFiles = inputSeedFile.read(); // don't pass pdb prefix, so it returns
    //StructureIterator iterator(seedFiles.first, 1, &seedFiles.second);
    FuseCandidateFinder fuseHelper;
    StructureCache *structures = graph.getStructures();
    
    int i = 0;
    for (string path: seedFiles.first) {
        //while (iterator.hasNext()) {
        //vector<Structure*> seedVec = iterator.next(); // only one structure
        Structure *seed = structures->getStructure(path);
        vector<Structure*> seedVec = { seed };
        if (fuseHelper.findNearbySeeds(realPeptide, seedVec, -5.0).size() == 0)
            continue;
        //Structure *seed = seedVec[0];
        if ((i++) % 100 == 0)
            cout << "Structure " << i << endl;
        if (seed->chainSize() > 1)
            continue;
        cout << "Loaded structure " << seed->getName() << endl;
        
        // Score the seed
        high_resolution_clock::time_point startTime = high_resolution_clock::now();
        auto result = scorer.score(seed);
        high_resolution_clock::time_point endTime = high_resolution_clock::now();
        double time = duration_cast<seconds>( endTime - startTime ).count();
        cout << "time: " << time << endl;
        totalTime += time;
        numTimes += 1;
        
        outputSeedFile.write(seed->getName(), seed->getChain(0).getID(), { writePerResidueScores(result) });
    }
    
    cout << "Done! Average time per seed: " << totalTime << "/" << numTimes << " = " << totalTime / numTimes << endl;
    
    delete target;*/
    
    // 1/23/19: Get scores for poses
    /*Structure *target = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/ATarget.pdb");
    Structure *seed = new Structure("/home/ifs-users/venkats/dtermen/1DUY_decorate_sample/test_pose_4U6Y.pdb");
    rmsdParams rParams(1.2, 15, 1);
    contactParams cParams;
    SequenceCompatibilityScorer scorer(target, rParams, cParams, fasstDB, 2, 2, 0.4, 0.05, 0.25, 1, 8000, 0.7); // new StructureCompatibilityScorer(target, fParams, rParams, cParams, fasstDB);
    // Read contact counts
    string contactCountsPath = basePath + "seed_scores/contact_scores/";
    for (int i = 1; i <= 25; i++) {
        scorer.readContactCounts(contactCountsPath + "contact_scores_" + to_string(i) + ".csv");
    }
    
    cout << "Loaded scorer" << endl;
    high_resolution_clock::time_point startTime = high_resolution_clock::now();
    auto result = scorer.score(seed);
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    cout << "Scoring took " << duration_cast<seconds>( endTime - startTime ).count() << " seconds" << endl;
    cout << "Result: ";
    for (auto it: result) {
        cout << it.first->getName() << it.first->getResidueIndex() << " = " << it.second << ", " << endl;
    }
    cout << endl;
    delete target;
    delete seed;*/
    
    // 1/17/19: Build adjacency list for all seeds
    /*SeedGraph graph;
    for (int i = 0; i < 20; i++) {
        cout << "Loading candidates from batch " << i << endl;
        FuseCandidateFile file(basePath + "overlaps_res2/results" + to_string(i) + ".txt");
        graph.load(file, dataPath);
    }
    
    cout << "Built graph with " << graph.seedSize() << " seeds and " << graph.residueSize() << " residues" << endl;
    
    cout << "Finding neighborhoods" << endl;
    vector<SeedGraph> subgraphs = graph.subgraphs();
    int subgraphIndex = 0;
    for (SeedGraph g: subgraphs) {
        g.write(outPath + "subgraph_" + to_string(subgraphIndex++) + ".txt");
    }*/
    
    // This writes fuse candidates given a seed list file (1/29/19: updated for general overlap detection
    // Multithreaded: use task ID 1-10
    /*float rmsdCutoff = 0.5f;
    //outPath = basePath + "results_general_0.5A/";
    SeedListFile seedFile(filesPath);
    auto fileContents = seedFile.read(dataPath);
    vector<string> candidatePaths = fileContents.first;
    //vector<string> chainIDs = fileContents.second;
    FuseCandidateFinder fuser(2, general, rmsdCutoff, 10, mode - 1);
    fuser.writeFuseCandidates(candidatePaths, nullptr, outPath);*/
    
    // 1/15/19: Try fusing actual fragments
    /*FuseCandidateFile file("/home/ifs-users/venkats/dtermen/nearby_seeds_test/new_results_res2/results0.txt");
    FuseHelper helper;
    helper.writeFuseResultsFromFile(file, "/home/ifsdata/scratch/grigoryanlab/venkats/dtermen/fuse_res2_test/", dataPath);*/
    
    /*const string basePath = "/home/ifs-users/venkats/dtermen/sandbox/";
    string path1 = basePath + "data/nchains1_nres5_score0.193_0.000_rmsd0.170_term097069_chidsA_match0.pdb";
    string path2 = basePath + "data/nchains1_nres5_score2.308_0.000_rmsd0.357_term004090_chidsA_match32.pdb";

    Structure s(path1, "");
    Structure s2(path2, "");
    for (Residue *r : s.getChain(0).getResidues()) {
        cout << "Residue: " << r->getName() << endl;
    }
    for (Residue *r : s2.getChain(0).getResidues()) {
        cout << "Residue 2: " << r->getName() << endl;
    }
    
    vector<vector<Residue *>> residues(18);
    residues[0].push_back(&s.getResidue(0));
    residues[1].push_back(&s.getResidue(1));
    residues[2].push_back(&s.getResidue(2));
    residues[3].push_back(&s.getResidue(3));
    residues[4].push_back(&s.getResidue(4));
    residues[3].push_back(&s2.getResidue(0));
    residues[4].push_back(&s2.getResidue(1));
    residues[5].push_back(&s.getResidue(5));
    residues[6].push_back(&s.getResidue(6));
    residues[7].push_back(&s.getResidue(7));
    residues[8].push_back(&s.getResidue(8));
    residues[9].push_back(&s.getResidue(9));
    residues[10].push_back(&s2.getResidue(2));
    residues[11].push_back(&s2.getResidue(3));
    residues[12].push_back(&s2.getResidue(4));
    residues[13].push_back(&s2.getResidue(5));
    residues[14].push_back(&s2.getResidue(6));
    residues[15].push_back(&s2.getResidue(7));
    residues[16].push_back(&s2.getResidue(8));
    residues[17].push_back(&s2.getResidue(9));
    fusionTopology topology(residues);
    cout << topology.length() << endl;
    cout << topology.numMobileAtoms() << endl;
    cout << topology.numChains() << endl;
    fusionOutput output;
    Fuser myFuser;
    Structure fusedStruct = myFuser.fuse(topology, output);
    cout << "Output: " << endl;
    cout << output.getScore() << endl;
    fusedStruct.writePDB(basePath + "data/test_fused.pdb");*/
    /*topology.addFragment(s);
    vector<int> overlapIndices { 4, 5, 11, 12, 13, 14, 15, 16, 17, 18 };
    topology.addFragment(s2, overlapIndices);*/
    
    return 0;
}
