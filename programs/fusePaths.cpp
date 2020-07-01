//
//  fusePaths.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/16/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
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
  opts.setTitle("Fuses a set of specified paths.");
  opts.addOption("target", "Path to the target PDB structure file", true);
  opts.addOption("peptideChain", "Chain ID for the peptide chain in the target, if one exists - it will be removed", false);
  opts.addOption("seeds", "Path to a binary file containing seed structures", true);
  opts.addOption("seedGraph","Path to a seed graph that contains the overlaps between seeds in the binary file",true);
  opts.addOption("seedChain", "Chain ID for the seed structures (default is '0')", false);
  opts.addOption("out", "Path to a directory into which the fused seed path structures and scores will be written", true);
  opts.addOption("paths","Path to a text file where each line specifies a path. Format: seed_A:residue_i,seed_B:residue_j,etc...",false);
  opts.addOption("config", "The path to a configfile",true);
  opts.setOptions(argc, argv);
  
  string targetPath = opts.getString("target");
  string binaryFilePath = opts.getString("seeds");
  string seedGraphPath = opts.getString("seedGraph");
  string outputPath = opts.getString("out");
  string pathsPath = opts.getString("paths");
  string configFilePath = opts.getString("config");
  
  // Load target structure and remove native peptide
  Structure target(targetPath);
  
  if (opts.isGiven("peptideChain")) {
    Chain *peptide = target.getChainByID(opts.getString("peptideChain", "B"));
    if (peptide != nullptr)
      target.deleteChain(peptide);
  }
  
  cout << "Loading seeds" << endl;
  StructuresBinaryFile* seedFile = new StructuresBinaryFile(binaryFilePath);
  seedFile->scanFilePositions();
  
  cout << "Loading graph" << endl;
  StructureCache* cache = new StructureCache(seedFile);
  SeedGraph* seedG = new SeedGraph(false,cache);

  
  if (!MstSys::fileExists(outputPath)) {
    MstSys::cmkdir(outputPath);
  }
  
  //read paths file
  vector<string> path_specifiers = MstUtils::fileToArray(pathsPath);
  
  // Scorer
  FragmentParams fParams(2, true);
  rmsdParams rParams(1.2, 15, 1);
  contactParams cParams;
  StructureCompatibilityScorer scorer(&target, fParams, rParams, cParams, configFilePath);
  
  ofstream out(MstSystemExtension::join(outputPath, "fused_paths.csv"), ios::out);
  if (!out.is_open())
    MstUtils::error("Could not open file stream");
  // CSV header
  out << "name,path,path_len,designability,num_contacts,num_designable,corroboration" << endl;
  
  //fuse specified paths
  cout << "trying to fuse paths..." << endl;
  SeedGraphPathSampler sampler(&target,seedG);
  vector<PathResult> pathResults = sampler.fusePaths(path_specifiers);
  
  cout << "Score paths" << endl;
  for (int i = 0; i < pathResults.size(); i++) {
    PathResult& result = pathResults[i];
    
    //get name
    string name = "fused_path_" + to_string(i);
    string name_whole = "fused_path_and_context_" + to_string(i);
    cout << "Scoring " << name << endl;
    
    out << name << ",";
    for (Residue *res: result.getOriginalResidues()) {
      out << res->getStructure()->getName() << ":" << res->getResidueIndex() << ";";
    }
    out << "," << result.size() << ",";
    
    //get the whole structure (including target-aligned residues) and path only
    Structure& whole = result.getFusedStructure();
    Structure path;
    result.getFusedPathOnly(path);
    
    // Score the path
    mstreal totalScore;
    int numContacts;
    int numDesignable;
    scorer.score(&path, totalScore, numContacts, numDesignable);
    cout << "Score: " << totalScore << endl;
    out << totalScore << "," << numContacts << "," << numDesignable << "," << corroborationScore(result, 2) << endl;
    
    // write out the structures
    path.setName(name);
    path.writePDB(MstSystemExtension::join(outputPath,name+".pdb"));
    
    whole.setName(name_whole);
    whole.writePDB(MstSystemExtension::join(outputPath,name_whole+".pdb"));
  }
  
  cout << "Done" << endl;
  out.close();
  
  return 0;
}

