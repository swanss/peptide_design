#ifndef coverage_h
#define coverage_h

#include "msttypes.h"
#include "mstsystem.h"
#include "mstcondeg.h"
#include "dtermen.h"
#include "mstmagic.h"
#include "mstfasst.h"

#include "pathsampler.h"
#include "structure_iter.h"
#include "termextension.h"
#include "utilities.h"
#include "mstsystem_exts.h"

using namespace MST;

class benchmarkUtils;

/* --------- seedSubstructureInfo --------- */

struct seedSubstructureInfo {
public:
    string structure_name = "";
    string chain_ID;
    // N-terminal residue within seed
    int res_idx;
    // number of residues in segment
    int res_length;
    //the rank of the match which the seed was generated from
    int match_number;
    //the RMSD of the protein segment(s) to the match
    mstreal match_rmsd;
    //indicates whether the match has the same amino acid at its central res
    bool sequence_match;
    mstreal rmsd;
    mstreal max_ca_deviation;
    mstreal alignment_cos_angle;
    
    seedSubstructureInfo() {};
    seedSubstructureInfo(string name, string ID, int idx, int length, int _match_number, mstreal _match_rmsd, bool _sequence_match, mstreal _rmsd, mstreal _alignment_cos_angle) : structure_name(name), chain_ID(ID), res_idx(idx), res_length(length), match_number(_match_number), match_rmsd(_match_rmsd), sequence_match(_sequence_match), rmsd(_rmsd), alignment_cos_angle(_alignment_cos_angle) {};
    
    bool operator < (const seedSubstructureInfo& other) const {
        return (rmsd < other.rmsd);
    }
    
    vector<Residue*> loadSeedSubstructureFromCache(StructureCache *seedCache, bool loadAnchor = false) {
        vector<Residue*> seedRes;
        string seedChainID = "0";
        Structure* seed = seedCache->getStructure(structure_name);
        Chain* seedChain = seed->getChainByID(seedChainID);
        for (int i = res_idx; i < res_idx+res_length; i++) {
            if (i < 0 || i >= seedChain->residueSize()) MstUtils::error("Seed residue ID not in seed chain","seedSubstructureInfo::loadSeedSubstructureFromCache");
            seedRes.push_back(&seedChain->getResidue(i));
        }
        if (loadAnchor) {
            for (int chain_idx = 0; chain_idx < seed->chainSize(); chain_idx++) {
                Chain* C = &seed->getChain(chain_idx);
                if (C->getID() == seedChainID) continue;
                vector<Residue*> anchorRes = C->getResidues();
                seedRes.insert(seedRes.end(),anchorRes.begin(),anchorRes.end());
            }
        }
        return seedRes;
    }
};

/* --------- sortedBins --------- */

class sortedBins {
public:
    sortedBins() {}
    
    sortedBins(mstreal _max_val, mstreal _interval = 0.1) {
        max_val = _max_val;
        interval = _interval;
        
        buildBins();
    }
    
    void insert(seedSubstructureInfo& entry, mstreal val);
    
    vector<seedSubstructureInfo> getBinByValue(mstreal val);
    vector<seedSubstructureInfo> getLowestValuePopulatedBin();
    bool allBinsEmpty() {return empty;}
    seedSubstructureInfo getLowestRMSDSeed(set<string> seedNames = {});
    //min_rmsd, max_rmsd, number of seeds aligned to segment
    vector<tuple<mstreal,mstreal,long>> getSeedsByBin();
    
    //  //min_rmsd, max_rmsd, number of seeds with a protein aligned region containing R_prot
    //  vector<tuple<mstreal,mstreal,long>> getNumFragmentsCoveringContact(int R_prot_idx);
    
    vector<seedSubstructureInfo> getAllSeeds();
    
    
    //sorts the entries within each bin by RMSD
    void sortBins();
    bool areBinsSorted() {return sorted;}
    
    void reset();
    
protected:
    void buildBins();
    int val2Bin(mstreal val) {
        return val / interval;
    };
    
    
    
private:
    vector<vector<seedSubstructureInfo>> bins;
    mstreal max_val, interval;
    bool sorted = true, empty = true;
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
    void setMaxMatchNumber(int _match_number_cutoff) {match_number_cutoff = _match_number_cutoff;}
    void setSeeds(string binFilePath);
    void setSeeds(StructuresBinaryFile* bin);
    
    //map seeds to peptide segments
    void findCoveringSeeds();
    bool mapSeedToChainSubsegments(vector<Atom*> seed_atoms, vector<Residue*> seed_residues, int match_number, mstreal match_rmsd, bool seq_match, bool only_check_if_aligned = false);
    
    //report the results
    //Info pertaining to the peptide (which is to be covered)
    void writePeptideResidues(string outDir);
    
    // Method rewritten to allow for writing any set of contacts
    void writeContacts(string outDir, set<pair<Residue*,Residue*>> provided_contact_residues = {}, map<string,contactList> provided_all_types_contacts = map<string,contactList>());
    
    //Info pertaining to coverage
    void writeAllAlignedSeedsInfo(string outDir);
    void writeBestAlignedSeeds(string outDir, int numSeeds, bool write_structures = false);
    void writeSegmentCoverage(string outDir);
    
    Structure* getTargetStructure() {return target;}
    vector<Residue*> getBindingSiteRes() {return bindingSiteRes;}
    set<pair<Residue*,Residue*>> getContactingResidues() {return contact_residues;}
    
    /**
    The binding site is defined as all residues with atoms within distanceCutoff (Å) from atoms of the peptide
    
    @param distanceCutoff the distance cutoff in angstroms used to define binding site residues
    @return The set of all binding site resides in S
    */
    set<Residue*> getBindingSiteResByDistance(mstreal distanceCutoff = 4.0);
    
    Chain* getPeptideChain() {return peptide_chain;}
    
    void resetBins();
    
    map<string,contactList> defineContacts(Structure& complex, vector<Residue*> peptide_residues);
    
protected:
    //methods called by constructor
    void setParams(string RL_path);
    
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
    
    mstreal getMaxRMSDForSegLength(mstreal max_rmsd, int segment_length) {
        return max_rmsd + 0.25*segment_length;
    }
    
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
    int match_number_cutoff;
    
    // peptide-protein complex
    Structure* original_complex;
    Structure complex; //stripped of sidechains
    Chain* peptide_chain;
    int peptide_first_residue_index_in_structure;
    
    // protein only structure
    Structure* target;
    vector<Chain*> target_chains; //pointers from complex
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
    string binFilePath = "";
    StructuresBinaryFile* seeds = nullptr;
    
    friend class pathFromCoveringSeeds;
};


class pathFromCoveringSeeds {
public:
    typedef vector<Residue*> residueSegment;
    
    pathFromCoveringSeeds(interfaceCoverage *_coverage, int _segmentLength = 3, bool _forceChimera = false) : segmentLength(_segmentLength), forceChimera(_forceChimera) {
        coverage = _coverage;
        seedCache = new StructureCache(coverage->seeds);
        terminalRes = segmentLength / 2;
        // get the peptide residue segments that will be covered
        vector<residueSegment> peptideResidueSegmentsVector = coverage->chainResidueSubsegments[segmentLength-1];
        peptideResidueSegments = set<residueSegment>(peptideResidueSegmentsVector.begin(),peptideResidueSegmentsVector.end());
    };
    
    ~pathFromCoveringSeeds() {
        delete seedCache;
    }
    
    vector<string> getCoveringPath(int maxSeedLength = 5, int minSeedLength = 3);
        
    map<int,vector<Residue*>> getCoveringResiduePath() {return coveringResiduesPath;}
    
    void writeCoveringSeeds(string outputPath);
protected:
    /**
     Finds the seed segment that 1)  covers the most of the not-yet-covered peptide segments and 2) has
     the lowest RMSD and returns the seed info as well as all of the covered segments
     
     Note: function ensures that the seed segment is from a seed that has not already been used to
     cover the peptide
     
     int terminalRes; this is equivalent to overlapLength in pathSampler
     // seed segments may be no longer than maxSeedLength
     // and no shorter than minSeedLength
     int maxSeedLength, minSeedLength;
     
     @param maxSeedLength the maximum length seed that will be used to
     @return a pair where the first element is all the info about the seed and the second
     element is the peptide residues that were covered by this seed
     */
    pair<seedSubstructureInfo,residueSegment> getBestCoveringSeed(int maxSeedLength, int minSeedLength);
    
    /**
     Given the vector of all peptide residue pointers that are covered by the seed segment, finds
     how many of the not-yet-covered peptide residues would be covered by this segment.
     
     @param peptideResiduesCoveredBySeedSegment a window of peptide residue pointers covered by the seed segment
     @return all residues that would be newly covered by the addition of this seed
     */
    
    set<residueSegment> getResidueSegmentsCoveredBySeed(vector<Residue*> peptideResiduesCoveredBySeedSegment, bool removeSeg = false);
    
    /**
     Selects seed residues from seeds in coveringSeeds and arranges them from the N to C terminus
     of the peptide that they cover, such that they can be properly fused by PathSampler
     
     @return a vector of residues that can converted to a path string
     */
    map<int,vector<Residue*>> getPathResidues();

private:
    interfaceCoverage *coverage;
    
    int terminalRes;
    
    /*
     The peptide residue segments are the elements that are covered by aligned seeds. The length
     of these segments is dictated by segmentLengths.
     */
    int segmentLength = 3;
    set<residueSegment> peptideResidueSegments;
    
    bool forceChimera; //if true, requires that every seed segment comes from a different seed
    
    set<string> seedNames;
    map<int,seedSubstructureInfo> coveringSeeds; //the key is the idx in the peptide chain
    map<int,vector<Residue*>> coveringResiduesPath; //pointers to residues in the seeds
    vector<string> coveringPathString;
    
    bool verbose = true;
    
    StructureCache* seedCache;
};

// Functions for the coverage benchmark
class coverageBenchmarkUtils {
public:
    /**
     Writes out a tsv file with the paired positions between the fused backbone and native peptide and rmsd per position. Also computes RMSD
     over the entire structure and reports this value.
     
     Note: assumes that even if the two structures have unequal numbers of residues, that they have complete backbones and the residue number
     can be used to unambiguously pair the residues.
     */
    mstreal static writeRMSDtoFile(string outputPath, Structure& fusedBackbone, Structure& nativePeptide);
    
    /**
     Writes out the contacts between all residues from fused chains (there could be more than one) in a similar format to the interfaceCoverage
     */
    void static writeContactstoFile(string outputPath, interfaceCoverage *IC, Structure& fusedPathandTarget, set<string> peptideChains, string rotLibFile);
    
    void static getResiduesFromMap(string resMapPath, Structure& structureA, Structure& structureB, vector<Residue*>& selectedResA, vector<Residue*>& selectedResB);
    
    /**
     Find the set optimally covering seeds and fuse together

     The set of covering seeds is obtained by the following algorithm
     1) Define the structural elements to be covered: in this case, overlapping 3-res windows (aka segments) of the peptide.
     2) Find the seed that covers the most windows. If there is a tie between two seeds, choose the seed with the lowest RMSD.
     3) If all segments are covered, terminate. Otherwise, return to step 2
     */
    void static fuseCoveringSeeds(interfaceCoverage* IC, bool force_chimera, int max_seed_length_fuse, string fusDir, string pdb_id, Structure& target, Structure& complex, bool two_step_fuse, string RL);
};

#endif
