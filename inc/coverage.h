//
//  compareterms.h
//  TPD_swans
//
//  Created by Sebastian Swanson on 1/29/19.
//

#ifndef coverage_h
#define coverage_h

#include "msttypes.h"
#include "mstsystem.h"
#include "mstcondeg.h"
#include "dtermen.h"
#include "mstmagic.h"
#include "mstfasst.h"

#include "utilities.h"
#include "termextension.h"
#include "structure_iter.h"

#include "Util.h"
#include "vdwRadii.h"

using namespace MST;

class benchmarkUtils;


/* --------- seedSubstructureInfo --------- */

struct seedSubstructureInfo {
public:
  string structure_name;
  string chain_ID;
  int res_idx;
  int res_length;
  mstreal rmsd;
  set<pair<int,int>> seed_protein_contacts;
  
  seedSubstructureInfo() {};
  seedSubstructureInfo(string name, string ID, int idx, int length, mstreal _rmsd, set<pair<int,int>> _seed_protein_contacts) : structure_name(name), chain_ID(ID), res_idx(idx), res_length(length), rmsd(_rmsd), seed_protein_contacts(_seed_protein_contacts) {};
  seedSubstructureInfo(const seedSubstructureInfo& other) : structure_name(other.structure_name), chain_ID(other.chain_ID), res_idx(other.res_idx), res_length(other.res_length), rmsd(other.rmsd), seed_protein_contacts(other.seed_protein_contacts) {};
  
  bool operator < (const seedSubstructureInfo& other) const {
    return (rmsd < other.rmsd);
  }
  
  string get_contacts_string() {
    stringstream ss;
    int count = 0;
    for (auto cont : seed_protein_contacts) {
      if (count != 0) ss << " ";
      ss << cont.first << "," << cont.second;
      count++;
    }
    ss << endl;
    return ss.str();
  };
};

/* --------- sortedBins --------- */

class sortedBins {
public:
  sortedBins() {};
  
  sortedBins(mstreal _max_val, mstreal _interval = 0.1) {
    max_val = _max_val;
    interval = _interval;
    
    buildBins();
  }
  
  void insert(seedSubstructureInfo& entry, mstreal val);
  
  vector<seedSubstructureInfo> getBinByValue(mstreal val);
  
  //min_rmsd, max_rmsd, number of seeds aligned to segment
  vector<tuple<mstreal,mstreal,long>> getSeedsByBin();
  
//  //min_rmsd, max_rmsd, number of seeds with a protein aligned region containing R_prot
//  vector<tuple<mstreal,mstreal,long>> getNumFragmentsCoveringContact(int R_prot_idx);
  
  vector<seedSubstructureInfo> getAllSeeds();
  
  
  //sorts the entries within each bin by RMSD
  void sortBins();
  
  void reset();
  
protected:
  void buildBins();
  int val2Bin(mstreal val) {
    return val / interval;
  };
  
  
  
private:
  vector<vector<seedSubstructureInfo>> bins;
  mstreal max_val, interval;
};

/* --------- allChainSubsegments --------- */

class allChainSegments {
public:
  allChainSegments() {};
  allChainSegments(Chain* peptide, Structure* target, int max_segment_length, mstreal max_rmsd, string s_cid, string rotLibPath);
  
  ~allChainSegments() {
    delete rotLib;
  };
  
  void mapSeedToChainSubsegments(vector<Atom*> seed_atoms);
  
  void resetBins();
  
  void writeSegmentCoverage(fstream& output);
//  void writeContactCoverage(fstream& output, set<pair<Residue*,Residue*>> contact_residues);
  
  void writeSeedsToFile(fstream& output);
  
protected:
  void mapSegmentToChainSubsegments(vector<Atom*> seed_segment, int seed_position, int length);
  
  set<pair<int,int>> getContacts(vector<Atom*> seed_segment);
  
private:
  vector<vector<AtomPointerVector>> chainSubsegments;
  vector<vector<sortedBins>> chainSubsegmentsBins;
  
  mstreal max_rmsd;
  int max_allowable_segment_length;
  RMSDCalculator rmsd_calc;
  
  //for contacts
  RotamerLibrary* rotLib;
  Structure* target; //only the protein atoms
  string s_cid;
  
};

/* --------- interfaceCoverage --------- */

class interfaceCoverage {
public:
  //constructors
  interfaceCoverage(Structure& S, string p_id, mstreal max_rmsd, int max_seed_length, string RL_path);
  
  //destructors
  ~interfaceCoverage();
  
  //methods
  void findCoveringSeeds(string _binFilePath);
  void writeCoverageToFiles(string outDir);
  
  void writeAllAlignedSeedstoFile(string outDir);
  
  Structure* getTargetStructure() {return target;}
  vector<Residue*> getBindingSiteRes() {return bindingSiteRes;}
  
protected:
  //methods called by constructor
  void setParams(mstreal rmsd, int max_segment_length, string RL_path);
  void defineCoverageElements();
  
  //this function assumes the peptide comes after the protein
  void prepareForTERMExtension();
  
//  set<int> getExtendedFragmentProteinResidueIdx(Structure* EF) {
//    set<int> protein_res;
//    for (int chain_idx = 0; chain_idx < EF->chainSize(); chain_idx++) {
//      Chain& C = (*EF)[chain_idx];
//      if (C.getID() != seed_chain_id) {
//        vector<Residue*> chain_res = C.getResidues();
//        for (Residue* R : chain_res) protein_res.insert(R->getNum());
//      }
//    }
//    return protein_res;
//  };
  
  vector<Atom*> getBackboneAtoms(Chain* C) {
    vector<Atom*> bb_atoms;
    vector<Residue*> chain_res = C->getResidues();
    for (Residue* R : chain_res) {
      vector<Atom*> R_atom = RotamerLibrary::getBackbone(R);
      bb_atoms.insert(bb_atoms.end(),R_atom.begin(),R_atom.end());
    }
    return bb_atoms;
  };
  
private:
  
  // global variables
  RotamerLibrary RL;
  string RL_path;
  
  // benchmark parameters
  mstreal cd_threshold, int_threshold, bbInt_cutoff;
  mstreal max_rmsd;
  int max_seed_length;
  string seed_chain_id;
  
  // peptide-protein complex
  Structure complex;
  Chain* peptide_chain;
  
  // protein only structure
  Structure* target;
  vector<Residue*> bindingSiteRes;
  
  // peptide-protein contacts
  map<string,contactList> all_types_contacts;

  // structural elements
  /* These are the indices of residues, residue-residue contacts, etc, at the peptide
   interface. They represent the basic units of structural information that could be useful
   for designing peptides. The goal of the coverage benchmark is to find fragments from
   other proteins that align to the interface i.e. "cover" the elements contained in the
   section that they align to.
   */
  vector<Residue*> peptide_residues;
  set<pair<Residue*,Residue*>> contact_residues; //the union of the contacts
  
  // peptide substructures
  allChainSegments peptideSubsegments;
  
  // seed structures binary file
  string binFilePath;
};

/* --------- benchmarkUtils --------- */


//class benchmarkUtils {
//public:
//  /* -- Functions for computing the similarity between vectors -- */
//  
//  /* Takes two amino acid distributions and computes the dot product between them. In practice, the
//   second distribution generally represents the "true distribution" and the first is some approximation
//   of it. */
//  static mstreal residueDistCosineSimilarity(vector<Residue*> approx_dist, vector<Residue*> true_dist);
//  
//  /* Takes two amino acid distributions and computes the Kullback-Leibler Divergence between them. In practice, the
//   second distribution generally represents the "true distribution" and the first is some approximation
//   of it. */
//  static mstreal residueDistKLD(vector<Residue*> approx_dist, vector<Residue*> true_dist);
//  
//};

///* Inherits from the dTERMen class to access methods for computing the energies conditioned on different
// structural properties. Each function can be used to score a residue from a seed, or a residue in a
// native binder that corresponds to a residue in a seed. */
//class seedScoring: public dTERMen {
//public:
//  /* default constructor */
//  seedScoring(const string& configFile);
//  
//  /* Calculates the probability distribution of amino acids at the given residue in the seed.
//   
//   In order to generate structural fragments to search against the PDB, the following is included:
//   -The central residue (cen_R)
//   -The residues +/- cen_res (specified by pmSeed)
//   -The residues contacting the central residue on the target protein (other chain)
//   
//   Much of this function is copied from dTERMen::buildEnergyTable
//   */
//  vector<mstreal> seedDist(Residue* cen_R);
//  
//  //    /* Includes: background, phi/psi, omega and own-backbone (self-residual). */
//  //    vector<mstreal> seedAloneDist(Residue* R);
//  
//  /* Includes: background */
//  vector<mstreal> backgroundDist(Residue* R);
//  
//  /* Includes: N/A. */
//  vector<mstreal> uniformDist();
//  
//  /* --- Getters --- */
//  vector<string> getTripleLetterAlph() {return triple_letter_alph;}
//  
//private:
//  FASST* F_p;
//  RotamerLibrary* RL_p;
//  vector<res_t> globalAlph;
//  vector<string> triple_letter_alph;
//  
//};

#endif
