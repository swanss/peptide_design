#ifndef Overlaps_H
#define Overlaps_H

#include <utility>
#include <set>
#include <map>
#include <unordered_map>

#include "msttypes.h"
#include "fileutilities.h"
#include "utilities.h"


// ============= Hasher ==============


/**
 An abstract class template for hashing residues and retrieving regions of
 hash values within a certain distance cutoff.
 */
template <typename H, typename D = MST::mstreal>
class ResidueHasher {
public:
    typedef H hash_type;
    typedef D distance_metric;
    
    /**
     Computes a hash value for the position and/or orientation of the given residue.
     
     @param res The residue to hash.
     
     @return A hash value in which to place the given residue.
     */
    virtual hash_type hash(MST::Residue *res) const;
    
    /**
     Determines a set of hash values whose buckets will comprise all residues within
     the given distance cutoff of the provided residue.
     
     @param res The residue to search around.
     @param cutoff The maximum distance around which to return hash values.
     
     @return A vector of hash values in which to find nearby residues.
     */
    virtual vector<hash_type> region(MST::Residue *res, distance_metric cutoff) const;
};

/**
 A concrete subclass of ResidueHasher that hashes residues by the position of
 the alpha carbon.
 
 Although this class is declared as a class template, it can only be instantiated
 with unsigned long as the hash type, and mstreal as the distance metric (see
 overlaps.cpp).
 */
template <typename H = unsigned long, typename D = MST::mstreal>
class CAResidueHasher: public ResidueHasher<H, D> {
public:
    typedef H hash_type;
    typedef D distance_metric;

    /**
     Initializes the hasher to hash positions within the given
     bounding box and with bins at the given increment.
     
     @param bbox A vector of six values, (xmin, xmax, ymin, ymax,
        zmin, zmax)
     @param increment The size of each bin in angstroms
     */
    CAResidueHasher(vector<MST::mstreal> bbox, MST::mstreal increment);
    
    hash_type hash(MST::Residue *res) const override;
    vector<hash_type> region(MST::Residue *res, distance_metric cutoff) const override;
    
private:
    vector<MST::mstreal> _bbox;
    MST::mstreal _increment;
    hash_type _xBinSize;
    hash_type _yBinSize;
    hash_type _zBinSize;
};


// ============= Verifier ==============


/**
 Abstract base class that exposes functionality to verify that a given set of
 residues do overlap.
 */
class OverlapVerifier {
public:
    virtual bool verify(const vector<MST::Residue *> &segment1, const vector<MST::Residue *> &segment2) const = 0;
    virtual ~OverlapVerifier() = 0;
};

/**
 Concrete subclass of OverlapVerifier that uses the maximum deviation between
 alpha carbons to check whether two segments are overlapping.
 */
class MaxDeviationVerifier: public OverlapVerifier {
public:
    MaxDeviationVerifier(MST::mstreal cutoff): _cutoff(cutoff) {};
    ~MaxDeviationVerifier() override {}
    
    bool verify(const vector<MST::Residue *> &segment1, const vector<MST::Residue *> &segment2) const override;

private:
    MST::mstreal _cutoff;
};

/**
 Concrete subclass of OverlapVerifier that uses the average cosine angle between
 normal vectors for each residue to check whether two segments are overlapping.
 */
class NormalVectorVerifier: public OverlapVerifier {
public:
    NormalVectorVerifier(MST::mstreal minCosAngle): _minCosAngle(minCosAngle) {};
    ~NormalVectorVerifier() override {}
    
    bool verify(const vector<MST::Residue *> &segment1, const vector<MST::Residue *> &segment2) const override;
    
private:
    MST::mstreal _minCosAngle;
};

/**
 Overlap verifier that combines the results of two child verifiers using an AND
 operation. Manages ownership of the child verifiers and deletes them upon destruction
 of this verifier.
 */
class CompositeVerifier: public OverlapVerifier {
public:
    CompositeVerifier(OverlapVerifier *v1, OverlapVerifier *v2): _v1(v1), _v2(v2) {};
    ~CompositeVerifier();
    
    bool verify(const vector<MST::Residue *> &segment1, const vector<MST::Residue *> &segment2) const override;

private:
    OverlapVerifier *_v1, *_v2;
};


// ============= Overlap Finder ==============


/**
 A helper data type that represents a position in a structure for possible overlap.
 */
struct OverlapCandidate {
    Structure *structure;
    string chainID;
    int offset;
    
    OverlapCandidate(Structure *s, string c, int o): structure(s), chainID(c), offset(o) {};
    
    /**
     Compares two overlap candidates to produce a deterministic ordering.
     
     @return true if the overlap candidates are in increasing order
     */
    static bool compare(const OverlapCandidate &s1, const OverlapCandidate &s2) {
        if (s1.structure != s2.structure)
            return s1.structure < s2.structure;
        if (s1.chainID != s2.chainID)
            return s1.chainID.compare(s2.chainID) < 0;
        return s1.offset < s2.offset;
    }
    
    /**
     @return true if the two candidates come from the same structure and chain
     */
    bool fromSameChain(const OverlapCandidate &other) const {
        return structure == other.structure && chainID == other.chainID;
    }
    
    /**
     Compares two overlap candidates for equality.
     */
    bool operator== (const OverlapCandidate &other) const {
        return (structure == other.structure &&
                chainID == other.chainID &&
                offset == other.offset);
    }

    bool operator!= (const OverlapCandidate &other) const { return !(*this == other); }
};

/**
 Class that uses residue hashing to efficiently find overlaps between seed structures.
 
 OverlapFinder utilizes two helper classes to modularize its functionality: the Hasher,
 which should inherit from ResidueHasher, and the OverlapVerifier. These classes permit
 the substitution of different overlap criteria, as long as the criterion can be expressed
 as a logical AND over all residues in the segment. For example:
 
 * The maximum distance between analogous C-alphas would be an appropriate criterion, since
   this is equivalent to requiring that the distance between C-alphas be at most some
   threshold value over all residues in the segment.
 * RMSD is NOT an appropriate criterion, since it cannot be linearly decomposed into each
   residue in a segment independently.
 
 In order to preserve the correctness of the algorithm, the Hasher's behavior must guarantee
 that for all pairs of segments that *would* be verified by the verifier, each residue of
 one segment appears within the `region` of the corresponding residue in the other segment.
 The region may of course return results that do not correspond to real overlaps, as these
 will be filtered out by the OverlapVerifier.
 */
template <class Hasher>
class OverlapFinder {
public:
    /**
     Initialize an OverlapFinder.
     
     @param hasher An instance of ResidueHasher that can hash and find regions using residues
     @param cutoff The distance cutoff to use for finding overlaps (passed to the ResidueHasher's
        region method)
     @param segmentLength The number of residues to use for overlap segments
     @param seedChain The chain ID corresponding to the seed (if empty, assumes the structures
        contain only the seeds)
     @param verifier An instance of OverlapVerifier that returns true if two sets of residues
        correspond to an overlap. If null, does not perform verification before returning
        overlaps.
     */
    OverlapFinder(const Hasher& hasher, typename Hasher::distance_metric cutoff, int segmentLength = 3, string seedChain = "", OverlapVerifier *verifier = nullptr): _hasher(hasher), _cutoff(cutoff), _segmentLength(segmentLength), _seedChain(seedChain), _verifier(verifier) {};
    
    /**
     Inserts the given set of structures into the hash map.
     */
    void insertStructures(const vector<MST::Structure *> &structures);
    
    /**
     Tests for overlaps between the given list of structures and the
     structures already present in the hash map.
     
     @param testStructures the structures to use as the "first structures" in
        an overlap
     @param constrainOrder Pass true if the testStructures and originally-hashed
        structures are the same set. This will avoid unnecessarily doubling
        returned overlaps
     */
    vector<FuseCandidate> findOverlaps(const vector<MST::Structure *> &testStructures, bool constrainOrder = false) const;
    
    /**
     Tests for overlaps using the same logic as the above call signature, but writes the
     results directly to a file.

     @param testStructures the structures to use as the "first structures" in
        an overlap
     @param outFile a file into which to write the results
     @param constrainOrder Pass true if the testStructures and originally-hashed
        structures are the same set. This will avoid unnecessarily doubling
        returned overlaps
     */
    void findOverlaps(const vector<MST::Structure *> &testStructures, FuseCandidateFile &outFile, bool constrainOrder = false) const;
    
    /**
     Finds overlaps with the given source of residues, and inserts them into the
     given results vector.
     */
    void findOverlaps(vector<MST::Residue *> &source, vector<FuseCandidate> &results, bool constrainOrder = false) const;
    
    /**
     Clears the current hash map contents.
     */
    void clear();
    
    // Indicates whether to print progress to the console.
    bool verbose = true;
    
private:
    int _segmentLength;
    typename Hasher::distance_metric _cutoff;
    string _seedChain;
    const Hasher &_hasher;
    OverlapVerifier *_verifier;
    
    // Maps hash values to vectors of pairs containing each residue and its index in the chain
    unordered_map<typename Hasher::hash_type, vector<pair<MST::Residue *, int>>> _hashMap;
        
    /**
     Determines whether the given overlap candidate is a real overlap, using the
     verifier if provided. If the overlap finder has no verifier, returns true.
     */
    bool checkOverlap(vector<MST::Residue *>::iterator begin, vector<MST::Residue *>::iterator end, const OverlapCandidate &candidate) const;
};

// This includes the template method implementations for OverlapFinder
#include "overlaps.tpp"

#endif
