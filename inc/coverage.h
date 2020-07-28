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
    mstreal alignment_cos_angle;
    
    seedSubstructureInfo() {};
    seedSubstructureInfo(string name, string ID, int idx, int length, mstreal _rmsd, mstreal _alignment_cos_angle) : structure_name(name), chain_ID(ID), res_idx(idx), res_length(length), rmsd(_rmsd), alignment_cos_angle(_alignment_cos_angle) {};
//    seedSubstructureInfo(const seedSubstructureInfo& other) : structure_name(other.structure_name), chain_ID(other.chain_ID), res_idx(other.res_idx), res_length(other.res_length), rmsd(other.rmsd) {};
    
    bool operator < (const seedSubstructureInfo& other) const {
        return (rmsd < other.rmsd);
    }
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
    vector<seedSubstructureInfo> getLowestValuePopulatedBin();
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

/* --------- interfaceCoverage --------- */

class interfaceCoverage {
public:
    //constructors
    interfaceCoverage(Structure& S, string p_id, string RL_path);
    
    //destructors
    ~interfaceCoverage();
    
    //set parameters
    void setSeqConst(bool val) {seq_const = val;}
    void setMatchRMSDConst(mstreal rmsd) {match_rmsd_cutoff = rmsd;}
    void setMaxRMSD(mstreal _max_rmsd) {max_rmsd = _max_rmsd;}
    void setMaxSegmentLength(int max_length) {max_allowable_segment_length = min(peptide_chain->residueSize(),max_length);}
    void setSeeds(string binFilePath);
    
    //map seeds to peptide segments
    void findCoveringSeeds();
    void mapSeedToChainSubsegments(vector<Atom*> seed_atoms, vector<Residue*> seed_residues);
    
    //report the results
    //Info pertaining to the peptide (which is to be covered)
    void writePeptideResidues(string outDir);
    void writeContacts(string outDir);
    
    //Info pertaining to coverage
    void writeAllAlignedSeedsInfo(string outDir);
    void writeBestAlignedSeeds(string outDir, int numSeeds);
    void writeSegmentCoverage(string outDir);
    
    Structure* getTargetStructure() {return target;}
    vector<Residue*> getBindingSiteRes() {return bindingSiteRes;}
    set<pair<Residue*,Residue*>> getContactingResidues() {return contact_residues;}
    
    void resetBins();
    
protected:
    //methods called by constructor
    void setParams(string RL_path);
    void defineCoverageElements();
    
    //this function assumes the peptide comes after the protein
    void prepareForTERMExtension();
    
    vector<Atom*> getBackboneAtoms(Chain* C) {
        vector<Atom*> bb_atoms;
        vector<Residue*> chain_res = C->getResidues();
        for (Residue* R : chain_res) {
            vector<Atom*> R_atom = RotamerLibrary::getBackbone(R);
            bb_atoms.insert(bb_atoms.end(),R_atom.begin(),R_atom.end());
        }
        return bb_atoms;
    };
    
    set<pair<int,int>> getContacts(vector<Atom*> seed_segment, int peptide_position);
    
    Structure* getSeedSegment(seedSubstructureInfo);
    
private:
    // global variables
    RotamerLibrary RL;
    string RL_path;
    
    // benchmark parameters
    mstreal cd_threshold, int_threshold, bbInt_cutoff;
    mstreal max_rmsd;
    int max_allowable_segment_length, max_seed_length;
    string seed_chain_id;
    bool seq_const;
    mstreal match_rmsd_cutoff;
    
    // peptide-protein complex
    Structure complex;
    Chain* peptide_chain;
    int peptide_first_residue_index_in_structure;
    
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
    
    // peptide segments that the seeds are mapped to
    vector<vector<vector<Residue*>>> chainResidueSubsegments;
    vector<vector<AtomPointerVector>> chainSubsegments;
    vector<vector<sortedBins>> chainSubsegmentsBins;
    
    RMSDCalculator rmsd_calc;
    
    // seed structures binary file
    string binFilePath;
    StructuresBinaryFile* seeds;
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
