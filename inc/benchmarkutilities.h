//
//  seedutilities.hpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 4/30/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include "msttypes.h"
#include "mstsystem.h"
#include "mstcondeg.h"
#include "dtermen.h"
#include "mstmagic.h"
#include "mstfasst.h"

#include "utilities.h"
#include "structure_iter.h"
#include "termextension.h"

class seedStatistics;

/* --------- histogram --------- */
class histogram {
public:
    histogram();
    histogram(mstreal _min_value, mstreal _max_value, int num_bins) : min_value(_min_value), max_value(_max_value) {
        bin_size = (max_value - min_value) / num_bins;
        bins.resize(num_bins);
    }
    histogram(string hist_file) {readHistFile(hist_file);}
    
    mstreal getMaxVal() {return max_value;}
    mstreal getMinVal() {return min_value;}
    mstreal getBinSize() {return bin_size;}
    mstreal getNumBins() {return bins.size();}
    
    //get the bin corresponding to the passed value and return its height
    mstreal getVal(mstreal value);
    mstreal getVal(int bin_id);
    
    void setBinVal(int bin_id, mstreal value) {
        if ((value < 0) || (value > bins.size())) MstUtils::error("Passed bin index outside of range for histogram: "+MstUtils::toString(bin_id)+" with max: "+MstUtils::toString(bins.size()));
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
    
private:
    vector<mstreal> bins;
    mstreal min_value;
    mstreal max_value;
    mstreal bin_size;
    
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
    rejectionSampler(histogram _proposal_hist, histogram _target_hist) : proposal_hist(_proposal_hist), target_hist(_target_hist) {
        if (target_hist.getNumBins() != proposal_hist.getNumBins()) MstUtils::error("Provided proposal and target histograms have a different number of bins");
        if (target_hist.getMinVal() != proposal_hist.getMinVal()) MstUtils::error("Provided proposal and target histograms have a different minimum values");
        if (target_hist.getBinSize() != proposal_hist.getBinSize()) MstUtils::error("Provided proposal and target histograms have a different bin widths");
        
        /* Checks if the values in the proposal histogram are always greater than the target histogram.
         If not, defines a constant, M, which scales the proposal histogram up. This ensures that the
         acceptance probability (target_prob / (M * proposal_prob)) is at most 1.
         */
        M = 1.0;
        for (int i = 0; i < proposal_hist.getNumBins(); i++) {
            if (target_hist.getVal(i) < min_value) target_hist.setBinVal(i,0.0);
            mstreal fraction = target_hist.getVal(i)/proposal_hist.getVal(i);
            if (fraction > M) M = fraction;
            if (M == std::numeric_limits<double>::max()) {
                string bin_num(MstUtils::toString(i));
                MstUtils::error("Proposal distribution scaling factor set to inf at bin  "+bin_num,"rejectionSampler::rejectionSampler");
            }
        }
        cout << "scaling the proposal distribution up by a factor of " << M << endl;
    };
    
    bool accept(mstreal value) {
        mstreal accept_prob = target_hist.getVal(value) / (proposal_hist.getVal(value) * M);
        mstreal sampled_value = MstUtils::randUnit();
        if (sampled_value <= accept_prob) return true;
        else return false;
    }
protected:
    mstreal getVal(mstreal value);
private:
    histogram proposal_hist;
    histogram target_hist;
    mstreal M; //the constant by which the proposal histogram should be multiplied by so it "envelopes" the target hist
    mstreal min_value = 0.000001; //any bin with a value below this in the target distribution is set to 0
};


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

/* --------- seedStatistics --------- */

class seedStatistics {
public:
    seedStatistics(Structure& S, string p_id, string seedBinaryPath = "");
    
    void setBinaryFile(string seedBinaryPath) {
        delete bin_file;
        bin_file = new StructuresBinaryFile(seedBinaryPath);
    }
    
    void writeStatisticstoFile(string output_path, string output_name, int num_final_seeds);
    
    histogram generateDistanceHistogram(mstreal min_value = 0, mstreal max_value = 25, int num_bins = 100, int sampled_seeds = 1000000);
    
    mstreal boundingSphereRadius(Structure* seed);
    mstreal centroid2NearestProteinAtom(Structure* seed);
    mstreal point2NearestProteinAtom(CartesianPoint point);
    //  mstreal atom2NearestProteinAtom(Structure* seed);
    
private:
    Structure& complex;
    Structure target;
    Chain* peptide;
    
    StructuresBinaryFile* bin_file;
    
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
    naiveSeedsFromBin(Structure& S, string p_id, string seedBinaryPath_in, string rotLibPath, rejectionSampler* sampler = nullptr);
    
    ~naiveSeedsFromBin() {
        delete target_PS;
    }
    
    /*
     Loads seeds from an existing seedBinaryFile and randomizes their position/orientation. During
     repositioning the seeds are placed within a bounding volume determined based on the original seed
     centroid distribution. If the seed placement results in clash, a new position/orientation are sampled,
     until there is no clash.
     */
    void newPose(string output_path, string out_name, bool position, bool orientation, vector<Residue*> binding_site = {});
    
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
    
    StructuresBinaryFile seeds;
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
    
    int max_attempts;
    mstreal distance;
    int neighbors;
    bool clash_check;
    bool rejection_sample;
};

/* --------- naiveSeedsfromDB --------- */

class naiveSeedsFromDB : public naiveSeedsFromBin {
public:
    naiveSeedsFromDB(Structure& S, string p_id, string seedBinaryPath_in, const string& dbFile, string rotLibPath, rejectionSampler* sampler = nullptr, int max_len = 50) : naiveSeedsFromBin(S,p_id,seedBinaryPath_in,rotLibPath,sampler), seedSampler(dbFile,max_len) {};
    
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

