#ifndef FindFuseable_H
#define FindFuseable_H

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <set>
#include <map>
#include <unordered_map>

#include "msttypes.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "seedutils.h"

using namespace std;
using namespace MST;

enum FuseCandidateSearchMode {
    /** searches all possible locations of overlap */
    general,
    /** searches overlap between C-terminus of seed and N-terminus of candidate */
    fromCTerm,
    /** searches overlap between N-terminus of seed and C-terminus of candidate */
    fromNTerm
};


class FuseCandidateFinder {
public:
    
    /** Initialize a FuseCandidateFinder with the given options */
    FuseCandidateFinder(int overlap = 1, FuseCandidateSearchMode mode = fromCTerm, float rmsd = 2.0f, int batchModulus = 1, int batchModValue = 0): numResOverlap(overlap), searchMode(mode), rmsdCutoff(rmsd), batchModulus(batchModulus), batchModValue(batchModValue) {}
    
    /**
     Identify seeds in the given list of structures whose bounding boxes
     intersect that of `seed`. Each Structure object is expected to contain only
     the chains in which overlap may occur - the `name` property of the Structure
     may be used to look up the full structure later.
     
     @param seed the seed to compare against
     @param structures the structures to filter
     @param inset an amount by which to inset the bounding box of seed (negative
            values increase the size of the bounding box)
     @return a vector of Structures whose bounding boxes overlap that of `seed`.
     */
    vector<Structure *> findNearbySeeds(Structure *seed, vector<Structure *> &structures, float inset = 0.0);

    /**
     Identify seeds in the given list of candidates that contain overlaps with
     `seed`. A candidate seed overlaps `seed` if for one of `seed`'s seed chains,
     the C terminus of `seed` has `numResOverlap` residues that align with the
     first `numResOverlap` residues of some seed chain in the candidate.
     
     Each Structure object is expected to contain only the chains in which
     overlap may occur. (The `name` property of the Structure may be used to
     look up the full structure later.)
     
     This method first filters for bounding box overlap using findNearbySeeds.
     
     @param seed the seed to compare against
     @param structures the structures to test, containing the seed chains only
     @param filterByBBox flag indicating whether to check all structures' bounding boxes before
        trying overlaps
     @return a vector of FuseCandidates that satisfy the overlap conditions above
     */
    vector<FuseCandidate> findOverlappingSeeds(Structure *seed, vector<Structure *> &structures, bool filterByBBox = true);
    
    /**
     * Finds overlaps using a proximity search object.
     */
    vector<FuseCandidate> findOverlappingSeeds(Structure *seed, vector<Structure *> &structures, ProximitySearch &ps);

    /**
     * Helper function that clients can use to quickly compute whether an RMSD is
     * less than the given threshold. Returns a large number if the RMSD is known
     * to be greater than the threshold.
     */
    mstreal quickRMSD(const vector<Atom *> &atoms1, const vector<Atom *> &atoms2, mstreal threshold);

    /**
     * Helper function that clients can use to check if a set of atoms forms an
     * overlap as defined by this fuse candidate finder. If the atom vectors are not
     * both of the expected length (numResOverlap * 4) or do not overlap, returns
     * a very large number. Overlap is defined by having an RMSD of at most rmsdCutoff
     * as well as an alignment at least minCosAngle (if provided).
     *
     * @param atoms1 the first set of atoms to check
     * @param atoms2 the second set of atoms to compare to
     * @param threshold the RMSD threshold to use. Uses the smaller of rmsdCutoff and this threshold.
     * @return the overlap RMSD if the atom vectors are the correct length and overlap, otherwise
     *  a very large number
     */
    mstreal atomsFormOverlap(const vector<Atom *> &atoms1, const vector<Atom *> &atoms2, mstreal threshold = -1.0);

    /**
     Compares all pairs of length-k residues from seedChain and candChain, where
     k = this->numResOverlap, and finds the minimum-RMSD overlap.
     
     @param baseChain the base (seed) chain
     @param candChain the candidate chain
     @return the best overlap starts for seedChain and candChain, and the RMSD
             associated with this comparison
     */
    vector<tuple<int, int, float>> scanOverlaps(Chain &baseChain, Chain &candChain);
    
    /**
     Checks that the Ca - Ca vectors for each set of atoms are aligned, i.e.
     their cos angle is > minCosAngle.
     
     @param atoms1 the set of atoms corresponding to the first chain
     @param atoms2 the set of atoms corresponding to the second chain
     @return true if the two chains go in roughly the same direction, and false otherwise
     */
    bool checkOverlapAlignment(const vector<Atom *> &atoms1, const vector<Atom *> &atoms2);
    
    /**
     Writes pairwise fuse candidates for every candidate in candidatePaths.
     
     @param candidatePaths the list of PDB file paths from which to identify
            pairwise candidates
     @param chainIDs the chain IDs to use as the seed chains within each PDB file (can be nullptr)
     @param outDir the directory in which to write the results, with trailing '/'
     */
    void writeFuseCandidates(vector<string> candidatePaths, vector<string> *chainIDs, string outDir, string baseSeedPath = "");
    /**
     Writes fuse candidates from the list of candidatePaths that could fuse
     with the seed at the given path.
     
     @param seedPath the PDB file to use as the seed from which to identify
            fuse candidates
     @param seedChainIDs the chain IDs to use within the above file
     @param candidatePaths the list of PDB file paths that could serve as the
            other part of a fuse candidate
     @param chainIDs the chain IDs to use for each PDB file in candidatePaths (can be nullptr)
     @param outPath the file path at which to write the results
     */
    void writeFuseCandidates(string seedPath, string seedChainIDs, vector<string> candidatePaths, vector<string> *chainIDs, string outPath, string baseSeedPath = "");
    
    /**
     Writes fuse candidates by reading structures from the given binary file.

     @param binaryFilePath a path to a .bin file containing seed structures
     @param outDir the directory at which to write the results
     */
    void writeFuseCandidates(string binaryFilePath, string outDir);

    /** the number of residues of the seed chains that must overlap */
    int numResOverlap = 1;
    /** the possible sites of overlap to use when generating fuse candidates */
    FuseCandidateSearchMode searchMode = fromCTerm;
    /** the maximum RMSD allowed between overlapping regions to generate a candidate */
    float rmsdCutoff = 2.0f;
    /** the minimum cos(angle) between Ca vectors allowed for overlaps */
    float minCosAngle = 0.0f;
    
private:
    bool intervalsOverlap(mstreal lo1, mstreal hi1, mstreal lo2, mstreal hi2);

    /** Writes pairwise fuse candidates to the given file. */
    void writeFuseCandidates(vector<Structure *> seeds, vector<Structure *> candidates, string outDir, string outFileName, string baseSeedPath = "", bool shouldAppend = false, int candidatesBatchIdx = -1);

    /** the mod to apply to the batch index when determining if the batch should be run, for writeFuseCandidates */
    int batchModulus;
    /** the value of the batch index mod batchModulus for the batch to be run, for writeFuseCandidates */
    int batchModValue;

    mstreal maxDeviation(const vector<Atom *> &s1, const vector<Atom *> &s2);
    double timeNew = 0.0;
    double timeOld = 0.0;
    int timeTrials = 0;

    unordered_map<int, vector<mstreal>> batchBoundingBoxes;
    vector<mstreal> computeBatchBoundingBox(vector<Structure *> seeds, double inset = 0.0);
};

#endif
