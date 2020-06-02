//
//  seedutilities.hpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 4/30/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#ifndef seedutilities_hpp
#define seedutilities_hpp

#include "msttypes.h"
#include "mstsystem.h"
#include "mstcondeg.h"
#include "dtermen.h"
#include "mstmagic.h"
#include "mstfasst.h"

#include "utilities.h"
#include "termextension.h"

#include "Util.h"
#include "vdwRadii.h"

/* --------- generateRandomSeed --------- */

class generateRandomSeed : public FASST {
public:
  /*
   Generates a seed structure by sampling a random window structures in the FASST DB.
   
   Note: the provided FASST DB must created with the option "s" in fasstdb.cpp, which splits
   chains with breaks into separate structures.
   */
  generateRandomSeed(const string& dbFile, int max_len);
  
  pair<Structure*,string> getSeed(int len);
private:
  int max_len;
  secondaryStructureClassifier classifier;
  
  /*
   In the map, the key is the length of the window to be sampled and the value is a vector with every
   possible window of that length.
   
   So, windows[l] is a vector, where every position specifies a unique window of length l.
   
   At each position of the vector there is a pair that specifies the target_index of the structure
   and the N_terminal residue index.
   
   windows[5][1000] is the thousandth length 5 window in the database.
   windows[5][1000].first is the target index of that window.
   windows[5][1000].second is the residue index of its N-terminus.
   
   windows are always contiguous sets of backbone atoms (as long as there are no breaks in the chain
   in the structures of the database).
   */
  map<int,vector<pair<int,int>>> windows;
};

/* --------- structureBoundingBox --------- */

struct structureBoundingBox {
public:
  //Use the just the backbone atoms to construct a bounding box, with the pad surrounding it
  structureBoundingBox(Chain* C, mstreal pad = 10.0);
  structureBoundingBox(vector<Residue*> residues, mstreal pad = 10.0);
  
  void construct_structureBoundingBox(AtomPointerVector atoms);
  
  mstreal xlo, xhi;
  mstreal ylo, yhi;
  mstreal zlo, zhi;
private:
  mstreal pad;
};

/* --------- naiveSeedsFromBin --------- */

class naiveSeedsFromBin {
  /*
   Generates a seed cloud that matches certain properties of existing clouds and saves to a new
   seedBinaryFile.
   */
public:
  naiveSeedsFromBin(Structure& S, string p_id, string seedBinaryPath_in, mstreal distance = 1.0, int neighbors = 1);
  
  ~naiveSeedsFromBin() {
    for (Atom* A : seed_centroids) delete A;
    delete target_PS;
    delete seed_PS;
  }
  
  /*
   Loads seeds from an existing seedBinaryFile and randomizes their position/orientation. During
   repositioning the seeds are placed within a bounding volume determined based on the original seed
   centroid distribution. If the seed placement results in clash, a new position/orientation are sampled,
   until there is no clash.
   */
  void newPose(string output_path, string out_name, bool position, bool orientation, vector<Residue*> binding_site = {});
  
protected:
  int transform(Structure* seed, structureBoundingBox& bounding_box, bool position, bool orientation, CartesianPoint new_centroid);
  
  Chain* peptide;
  StructuresBinaryFile seeds;
private:
  Structure& complex;
  Structure target;
  
  
  // for identifying seeds with clashes to the protein
  Structure target_BB_structure;
  AtomPointerVector target_BB_atoms;
  ProximitySearch* target_PS;
  
  // for identifying allowable centroid positions
  AtomPointerVector seed_centroids;
  ProximitySearch* seed_PS;
  
  int max_attempts;
  mstreal distance;
  int neighbors;
};

/* --------- naiveSeedsfromDB --------- */

class naiveSeedsFromDB : public naiveSeedsFromBin {
public:
  naiveSeedsFromDB(Structure& S, string p_id, const string& dbFile, string seedBinaryPath_in, mstreal distance = 1.0, int max_len = 50) : naiveSeedsFromBin(S,p_id,seedBinaryPath_in,distance), seedSampler(dbFile,max_len) {};
  
  /*
   Loads each seed from an existing seedBinaryFile, finds its residue length, and samples a new one
   from the DB and randomizes its position/orientation. During repositioning the seeds are placed
   within a bounding volume determined based on the original seed centroid distribution. If the seed
   placement results in clash, a new position/orientation are sampled, until there is no clash.
   */
  void newPose(string output_path, string out_name, bool position, bool orientation, vector<Residue*> binding_site = {});
  
  
private:
  generateRandomSeed seedSampler;
};

/* --------- seedStatistics --------- */

class seedStatistics {
public:
  seedStatistics(Structure& S, string p_id);
  
  void writeStatisticstoFile(string seedBinaryPath_in, string output_path, string output_name, int num_final_seeds);
  
  mstreal boundingSphereRadius(Structure* seed);
  mstreal centroid2NearestProteinAtom(Structure* seed);
//  mstreal atom2NearestProteinAtom(Structure* seed);
  
private:
  Structure& complex;
  Structure target;
  Chain* peptide;
  
  // variables stored for identifying seeds with clashes during randomization
  Structure target_BB_structure;
  AtomPointerVector target_BB_atoms;
  ProximitySearch* target_PS;
  
  mstreal neigborhood;
};

#endif /* seedutilities_hpp */
