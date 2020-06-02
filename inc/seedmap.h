//
//  seedmap.h
//  dtermen
//
//  Created by Venkatesh Sivaraman on 1/31/19.
//

#ifndef seedmap_h
#define seedmap_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <set>
#include <unordered_set>
#include <unordered_map>

#include "msttypes.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "seedutils.h"
#include "seedscore.h"

using namespace std;
using namespace MST;

/**
 This class manages a large list of seed structures in relation to a target
 structure, and allows for easy saving and retrieval of information about them
 such as contacts with the target, proximity to other seeds, and scoring.
 
 In this class, we refer to the "match" as the entire PDB structure of the protein
 from which a seed is derived. A match structure's getName() should uniquely
 identify it among other structures in the PDB, but other seeds may come from
 the same match.
 
 A "seed" is the PDB structure of the part aligned to the target, possibly
 including fragments directly overlapping with the target. This class maintains
 a mapping between seeds and their original match structures in order to share
 data for each residue in the match structure.
 */
class SeedMap {
public:
    SeedMap(Structure *target): _target(target) {};
    /**
     Initialize a seed map by reading it in from the given directory.
     @param directory the directory to read from (should have been written by
            another SeedMap instance)
     */
    SeedMap(string directory);
    
    /**
     Writes this seed map to a directory, creating the directory if necessary.
     @param directory the directory path to write to
     */
    //void write(string directory);
    
    /**
     Registers the given seed and its residues for saving.
     
     @param seed the seed structure (IMPORTANT: getName should be the destination
            file name of the seed)
     @param match the *entire* match structure from which seed is
            a subset
     @param resIndexes a vector where resIndexes[i] is the index in `match` of
            the analog to seed.getResidue(i), or -1 if seed.getResidue(i) is
            aligned to the target
     */
    //void addSeed(Structure *seed, Structure *match, vector<int> resIndexes);
    
    /**
     Registers the given score for the given contact in the seed. Also saves
     that a contact exists between the two residues, if not already saved. The
     seed must have been added before by a call to addSeed.
     
     @param matchRes the contacting residue from the PDB match
     @param targetRes the residue on the target
     @param score the score to save
     */
    //void addScore(Residue *matchRes, Residue *targetRes, mstreal score);
    
    /**
     Registers the given score for the given contact in the seed. Also saves
     that a contact exists between the two residues, if not already saved. The
     seed must have been added before by a call to addSeed.
     
     @param targetRes the residue on the target
     @param score the score to save
     */
    //void addTargetScore(Residue *targetRes, mstreal score);

    /**
     Registers that a contact exists between the two given residues.
     
     @param matchRes the residue on the seed (from the original PDB match)
     @param targetRes the residue on the target
     */
    //void addContact(Residue *matchRes, Residue *targetRes);
    
    /**
     Registers that the seed identified by seedId will be saved at savePath.
     
     @param residues the residues FROM THE MATCH STRUCTURE that are stored in
            this seed
     @param savePath the path to which the seed will be written
     */
    //void addSavePath(vector<Residue *> residues, string savePath);
    
    /// Sets the score types that will be saved.
    void setScoreIds(vector<string> scoreIds) { _scoreIds = scoreIds; }
    /// Gets the score types that will be saved.
    vector<string> getScoreIds() { return _scoreIds; }

    /// Sets the score types that will be saved for the target.
    void setTargetScoreIds(vector<string> scoreIds) { _targetScoreIds = scoreIds; }
    /// Gets the score types that will be saved.
    vector<string> getTargetScoreIds() { return _targetScoreIds; }

    
private:
    Structure *_target = nullptr;
    PositionMap<string> _residues;
    vector<string> _scoreIds;
    vector<string> _targetScoreIds;
    unordered_map<string, pair<string, int>> _seedOrigins;
    unordered_map<string, vector<string>> _contactsFromTarget;
    unordered_map<string, vector<string>> _contactsFromSeeds;
    unordered_map<string, unordered_map<string, mstreal>> _seedScores;
    unordered_map<string, unordered_map<string, mstreal>> _targetScores;

    void readMetadata(string inPath);
    void readTargetContacts(string inPath);
    void readSeedProximities(string inPath);
    void readResidueScores(string inPath);
    void readTargetScores(string inPath);
};
#endif /* seedmap_h */
