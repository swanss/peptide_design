#ifndef benchmarkutilities_h
#define benchmarkutilities_h

#include "msttypes.h"
#include "mstsystem.h"
#include "mstcondeg.h"
#include "dtermen.h"
#include "mstmagic.h"
#include "mstfasst.h"

#include "pathsampler.h"
#include "secondarystructure.h"
#include "structure_iter.h"
#include "utilities.h"
#include "vdwRadii.h"

class seedStatistics;

/* --------- histogram --------- */
class oneDimBinnedData {
public:
    oneDimBinnedData() {};
    oneDimBinnedData(mstreal _min_value, mstreal _max_value, int num_bins) : min_value(_min_value), max_value(_max_value) {
        bin_size = (max_value - min_value) / num_bins;
        bins.resize(num_bins);
    }
    oneDimBinnedData(string hist_file) {readHistFile(hist_file);}
    
    mstreal getMaxVal() {return max_value;}
    mstreal getMinVal() {return min_value;}
    mstreal getBinSize() {return bin_size;}
    mstreal getNumBins() {return bins.size();}
    
    //get the bin corresponding to the passed value and return its height
    mstreal getVal(mstreal value);
    mstreal getVal(int bin_id);
    
    void setBinVal(int bin_id, mstreal value) {
        if ((bin_id < 0) || (bin_id >= bins.size())) MstUtils::error("Passed bin index outside of range for histogram: "+MstUtils::toString(bin_id)+" with max: "+MstUtils::toString(bins.size()));
        bins[bin_id] = value;
    }
    
    /*
     The format is:
     lower_bound,upper_bound,value
     0,1,1
     1,2,3,
     ...
     */
    
    void readHistFile(string hist_file);
    
    void writeHistFile(string hist_file);
    
    void printInfo() {
        cout << "lower_bound,upper_bound,value" << endl;
        mstreal lower_bound = min_value;
        mstreal upper_bound = min_value + bin_size;
        for (mstreal value : bins) {
            cout << lower_bound << "," << upper_bound << "," << value << endl;
            lower_bound += bin_size;
            upper_bound += bin_size;
        }
    };
    
    void exceptOutOfRangeQueries() {outOfRangeException = true;}
    
private:
    vector<mstreal> bins;
    mstreal min_value;
    mstreal max_value;
    mstreal bin_size;
    
    bool outOfRangeException = false;
};

/* --------- rejectionSampler --------- */

class rejectionSampler {
    
    /*
     Rejection sampling allows us to sample from a non-parametric probability density (the target
     distribution) by first sampling from an easy-to-sample-from distribution (the proposal distribution)
     and then accept or rejecting, based on the relative density of the two distributions at the
     sampled value.
     
     For more info, check out:
     https://people.eecs.berkeley.edu/~jordan/courses/260-spring10/lectures/lecture17.pdf
     */
    
public:
    rejectionSampler(oneDimBinnedData _proposal_hist, oneDimBinnedData _target_hist) : proposal_hist(_proposal_hist), target_hist(_target_hist) {
        if (target_hist.getNumBins() != proposal_hist.getNumBins()) MstUtils::error("Provided proposal and target histograms have a different number of bins");
        if (target_hist.getMinVal() != proposal_hist.getMinVal()) MstUtils::error("Provided proposal and target histograms have a different minimum values");
        if (target_hist.getBinSize() != proposal_hist.getBinSize()) MstUtils::error("Provided proposal and target histograms have a different bin widths");
        
        /* Checks if the values in the proposal histogram are always greater than the target histogram.
         If not, defines a constant, M, which scales the proposal histogram up. This ensures that the
         acceptance probability (target_prob / (M * proposal_prob)) is at most 1.
         */
        M = 1.0;
        for (int i = 0; i < proposal_hist.getNumBins(); i++) {
            if (proposal_hist.getVal(i) == 0.0) continue; //don't divide by zero
            mstreal fraction = target_hist.getVal(i)/proposal_hist.getVal(i);
            if (fraction > M) M = fraction;
            if (isinf(M)) {
                string bin_num(MstUtils::toString(i));
                MstUtils::error("Proposal distribution scaling factor set to inf at bin  "+bin_num,"rejectionSampler::rejectionSampler");
            }
        }
        cout << "scaling the proposal distribution up by a factor of " << M << endl;
    };
    
    /*
     Note: if the value falls in a bin where the proposal distribution is set to zero, then it
     will always be rejected. Theoretically this is an issue because it means we cannot recapitulate
     the target distribution tails (bins with non-zero values will not be sampled), but in practice this
     shouldn't be a huge issue.
     */
    bool accept(mstreal value) {
        if (proposal_hist.getVal(value) == 0.0) return false;
        mstreal accept_prob = target_hist.getVal(value) / (proposal_hist.getVal(value) * M);
        mstreal sampled_value = MstUtils::randUnit();
        if (sampled_value <= accept_prob) return true;
        else return false;
    }
protected:
    mstreal getVal(mstreal value);
private:
    oneDimBinnedData proposal_hist;
    oneDimBinnedData target_hist;
    mstreal M; //the constant by which the proposal histogram should be multiplied by so it "envelopes" the target hist
};


/* --------- generateRandomSeed --------- */

class generateRandomSeed : public FASST {
public:
    /*
     Generates a seed structure by sampling a random window structures in the FASST DB.
     
     Note: the provided FASST DB must created with the option "s" in fasstdb.cpp, which splits
     chains with gaps. Otherwise, sampled seeds will have gaps.
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
     
     windows[5][1000] is the 1001th length 5 window in the database.
     windows[5][1000].first is the target index of that window.
     windows[5][1000].second is the residue index of its N-terminus.
     
     windows are always contiguous sets of backbone atoms (as long as there are no breaks in the chain
     in the structures of the database).
     */
    map<int,vector<pair<int,int>>> windows;
    
    /*
     num_sampled keeps track of how many times a window has been sampled. This is useful for ensuring
     that the sampled seeds always have a unique name, which is important for other applications.
     */
    map<int,vector<int>> num_sampled;
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

/* --------- seedStatistics --------- */

class seedStatistics {
public:
    seedStatistics(Structure& S, string p_id, string _seedBinaryPath);
    seedStatistics(Structure& S, string p_id, StructuresBinaryFile* _seed_bin = nullptr);
    
    ~seedStatistics() {
        if (seedBinaryPath != "") delete seed_bin;
    }
    
    void setBinaryFile(string _seedBinaryPath) {
        if (seedBinaryPath != "") delete seed_bin;
        seedBinaryPath = _seedBinaryPath;
        seed_bin = new StructuresBinaryFile(seedBinaryPath);
    }
    
    void setBinaryFile(StructuresBinaryFile* bin) {
        if (seedBinaryPath != "") delete seed_bin;
        seedBinaryPath = "";
        seed_bin = bin;
    }
    
    void writeStatisticstoFile(string output_path, string output_name, int num_final_seeds);
    
    oneDimBinnedData generateDistanceHistogram(mstreal min_value = 0, mstreal max_value = 50, int num_bins = 200, int sampled_seeds = 1000000);
    
    mstreal boundingSphereRadius(Structure* seed);
    mstreal centroid2NearestProteinAtom(Structure* seed);
    mstreal point2NearestProteinAtom(CartesianPoint point);
    //  mstreal atom2NearestProteinAtom(Structure* seed);
    
protected:
    void constructor(string p_id);
    
private:
    Structure& complex;
    Structure target;
    Chain* peptide;
    
    string seedBinaryPath = "";
    StructuresBinaryFile* seed_bin = nullptr;
    
    // variables stored for identifying seeds with clashes during randomization
    Structure target_BB_structure;
    AtomPointerVector target_BB_atoms;
    ProximitySearch* target_PS;
    
    mstreal neigborhood;
};

/* --------- naiveSeedsFromBin --------- */

class naiveSeedsFromBin {
    friend class naiveSeedsFromDB;
    /*
     Generates a seed cloud that matches certain properties of existing clouds and saves to a new
     seedBinaryFile.
     */
public:
    naiveSeedsFromBin(Structure& S, string p_id, StructuresBinaryFile* _bin, string rotLibPath, rejectionSampler* sampler = nullptr);
    
    ~naiveSeedsFromBin() {
        delete target_PS;
    }
    
    /*
     Loads seeds from an existing seedBinaryFile and randomizes their position/orientation. During
     repositioning the seeds are placed within a bounding volume determined based on the original seed
     centroid distribution. If the seed placement results in clash, a new position/orientation are sampled,
     until there is no clash.
     
     if num_seeds = 0, then will generate the same number of seeds in as input binary file
     */
    void newPose(string output_path, string out_name, bool position, bool orientation, int num_seeds = 0);
    
    void setRejectionSampler(rejectionSampler* _sampler) {sampler = _sampler;}
    
    void setClashChecking(bool _clash_check) {clash_check = _clash_check;}
    void setRejectionSampling(bool _rejection_sample) {
        rejection_sample = _rejection_sample;
        rejection_sample = true;
    }
    
protected:
    int transform(Structure* seed, structureBoundingBox& bounding_box, bool position, bool orientation, CartesianPoint new_centroid);
    
private:
    Structure& complex;
    Structure target;
    Chain* peptide;
    vector<Residue*> peptide_interface_residues;
    string seed_chain_id = "0";
    
    StructuresBinaryFile* seeds = nullptr;
    set<string> int_properties;
    set<string> real_properties;
    
    // for identifying seeds with clashes to the protein
    Structure target_BB_structure;
    AtomPointerVector target_BB_atoms;
    ProximitySearch* target_PS;
    
    //  // for identifying allowable centroid positions
    //  AtomPointerVector seed_centroids;
    //  ProximitySearch* seed_PS;
    
    // for sampling centroid positions
    rejectionSampler* sampler;
    seedStatistics stat;
    
    int max_attempts_ori;
    int max_attempts_pos;
    mstreal distance;
    int neighbors;
    bool clash_check;
    bool rejection_sample;
};

/* --------- naiveSeedsfromDB --------- */

class naiveSeedsFromDB : public naiveSeedsFromBin {
public:
    naiveSeedsFromDB(Structure& S, string p_id, StructuresBinaryFile* seed_bin, const string& dbFile, string rotLibPath, rejectionSampler* sampler = nullptr, int max_len = 50) : naiveSeedsFromBin(S,p_id,seed_bin,rotLibPath,sampler), seedSampler(dbFile,max_len) {};
    
    /*
     Loads each seed from an existing seedBinaryFile, finds its residue length, and samples a new one
     from the DB and randomizes its position/orientation. During repositioning the seeds are placed
     within a bounding volume determined based on the original seed centroid distribution. If the seed
     placement results in clash, a new position/orientation are sampled, until there is no clash.
     
     if num_seeds = 0, then will generate the same number of seeds in as input binary file
     */
    void newPose(string output_path, string out_name, bool position, bool orientation, int num_seeds = 0);
    
    
private:
    generateRandomSeed seedSampler;
};

/* --------- searchInterfaceFragments --------- */

/*
 This class is modeled after the coverage algorithm described in Vanhee et al. 2009
 */

struct interfaceFragment {
    Structure s;
    string name;
    pair<Residue*,Residue*> contact; //peptide,protein
    vector<int> fragResIdx;
    vector<string> cIDs; //the chain names of each chain in the fragment
    vector<Atom*> protein_atoms;
    vector<Atom*> peptide_atoms;
    int peptide_position;
    int segment_length;
    int matches;
    
    void reportFragment(string base_path) {
        cout<< "Name: " << name << endl;
        cout << "Structure has: " << s.chainSize() << " chains, " << s.residueSize() <<  " residues, and " << s.atomSize() << " atoms" << endl;
        cout << "Centered on contact between " << contact.first->getChainID() << contact.first->getNum() << " and " << contact.second->getChainID() << contact.second->getNum() << endl;
        cout << "Chain IDs: ";
        for (string chain_id : cIDs) cout << " " << chain_id;
        cout << endl;
        cout << "Protein atoms: " << protein_atoms.size() << endl;
        cout << "Peptide atoms: " << peptide_atoms.size() << endl;
        cout << "Writing fragment..." << endl;
        s.writePDB(base_path + name+".pdb");
    }
    //set structure name
    /*
     e.g. 1A1M_A_C_1_1-A59-A171A63A167-2_2-0.528887-1
     structure name
     central residue (protein)
     contacting residue (peptide)
     flanking residues (protein)
     flanking residues (peptide)
     */
    void setName(string complex_name, int flank) {
        name += complex_name + "-";
        name += contact.second->getChainID() + MstUtils::toString(contact.second->getNum()) + "-";
        name += contact.first->getChainID() + MstUtils::toString(contact.first->getNum()) + "-";
        name += MstUtils::toString(flank) + "_" + MstUtils::toString(flank);
    }
};

class searchInterfaceFragments {
public:
    searchInterfaceFragments(set<pair<Residue*,Residue*>> contact_residues, string fasstDBpath);
    
    void findMatches(string base_path);
    
    void setFlankingResidues(int _flank) {flank = _flank;}
    
    void writeFragments(string base_path);
    
protected:
    /*
     Checks if the given residues have enough flanking residues to create a proper fragment
     {5,5}. If not, the pointer is reassigned to a new residue that has sufficient flanking residues
     (but is still close enough to contain the passed residues in the final fragment). Next checks
     if the peptide/protein residue have already been added to the set, if not, adds and returns
     "true"
     */
    bool checkViableContact(Residue* R_pep, Residue* R_prot);
    
    string getNamePrefix(const interfaceFragment& f);
    
    void renameChains(Structure& match, const interfaceFragment& f);
    
    void repositionMatch(Structure& match, const interfaceFragment& f, mstreal& protein_rmsd_before_realign, mstreal& protein_rmsd_after_realign, mstreal& peptide_rmsd_before_realign, mstreal& peptide_rmsd_after_realign);
    
    vector<Atom*> getProteinAlignedAtoms(const Structure& match, const interfaceFragment& f, vector<Atom*>& peptide_aligned_atoms);
    
    /* Assumes that all fragments are two segment */
    mstreal segmentCentroidDistances(const Structure& match_split);
    
private:
    // Main variables and parameters
    Structure* complex;
    string complex_name;
    set<pair<Residue*,Residue*>> contact_residues; //peptide,protein
    set<pair<Residue*,Residue*>> interface_central_residues; //subset of contact_residues used to define fragments
    string peptide_cid;
    string seed_cid;
    
    vector<interfaceFragment> interface_fragments;
    
    // Fragment definition
    int flank;
    
    // Search parameters
    FASST F;
    fasstSearchOptions foptsBase; // base FASST options that will be used with every search
    string fasstdbPath;
    mstreal max_rmsd;
    int top_N_matches;

    MstTimer timer;
};

#endif
