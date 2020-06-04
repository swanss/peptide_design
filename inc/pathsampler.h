//
//  pathsampler.h
//  dtermen
//
//  Created by Venkatesh Sivaraman on 1/24/20.
//

#ifndef pathsampler_h 
#define pathsampler_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_set>
#include <unordered_map>

#include "msttypes.h"
#include "mstcondeg.h"
#include "mstfuser.h"
#include "mstsystem_exts.h"
#include "seedgraph.h"
#include "clustertree.h"
#include "Util.h" // For pair_hash

class PathResult {
public:
    vector<Residue *> getOriginalResidues() { return _originalResidues; }
    int size() { return _originalResidues.size(); }

    Structure& getFusedStructure() { return _fusedPath; }
    void getFusedPathOnly(Structure &ret) {
        Chain *c = new Chain;
        for (int i = _seedStartIdx; i < _fusedPath.residueSize(); i++) {
            c->appendResidue(new Residue(_fusedPath.getResidue(i)));
        }
        ret.appendChain(c);
    }
    PathResult(vector<Residue *> originalResidues, Structure fusedPath, int seedStartIdx): _originalResidues(originalResidues), _fusedPath(fusedPath), _seedStartIdx(seedStartIdx) {}

private:
    vector<Residue *> _originalResidues;
    Structure _fusedPath;
    int _seedStartIdx;
};

class PathSampler {
public:
    PathSampler(Structure *target, SeedGraph *graph, StructureCache *structures): _target(target), _graph(graph), _structures(structures), targetAPV(target->getAtoms()), ps(targetAPV, 10.0) {};

    PathSampler(Structure *target, FragmentFetcher *fetcher, ClusterTree *overlapTree, int overlapSize, mstreal overlapRMSD, vector<Residue *> residues): _target(target), _fetcher(fetcher), _searchTree(overlapTree), _searchTreeResidues(residues), overlapLength(overlapSize), overlapRMSD(overlapRMSD), targetAPV(target->getAtoms()), ps(targetAPV, 10.0) {};

    vector<PathResult> sample(int numPaths);
    
    // If true, constrain sampled paths to never use seeds that have been used in previously-sampled paths
    bool uniqueSeeds = false;

    // If non-null, weight possible extensions by the predominant secondary structure classification of the seed
    string *preferredSecondaryStructure = nullptr;
    float secondaryStructureWeight = 2.0f;

private:
    Structure *_target = nullptr;
    SeedGraph *_graph = nullptr;
    StructureCache *_structures = nullptr;
    FragmentFetcher *_fetcher = nullptr;
    ClusterTree *_searchTree = nullptr;
    vector<Residue *> _searchTreeResidues;

    // Used if uniqueSeeds is true
    unordered_set<string> _usedSeeds;

    // Number of residues to include for overlap when fusing
    int overlapLength = 3;

    mstreal overlapRMSD = 0.0;

    AtomPointerVector targetAPV;
    ProximitySearch ps;

    // Minimum cosine angle between first and last alpha carbons in overlap
    mstreal minCosineAngle = 0.75;

    void sampleFromGraph(int numPaths, vector<PathResult> &results);
    void sampleFromOverlapTree(int numPaths, vector<PathResult> &results);

    mstreal cosineAngle(const vector<Residue *> &res1, const vector<Residue *> &res2);

    vector<int> shuffleResultIndexes(ClusterSearchResults &searchResults);

    /**
     * Generates a PathResult by fusing the given path together, making 
     * sure it is valid, and checking for clashes.
     */
    bool emplacePathFromResidues(vector<Residue *> path, vector<PathResult> &results);

    AtomPointerVector buildAPV(const vector<Residue *>::const_iterator &begin, const vector<Residue *>::const_iterator &end);

    /**
     * Gets the residues from the target that correspond to the non-seed
     * match residues in the given path.
     *
     * @param residues a list of path residues
     * @return a list of residues from _target that are marked as associated
     *  with the seeds that the path residues come from
     */
    vector<Residue *> getTargetResidues(const vector<Residue *> &residues);
    /**
     * Gets the match residues (which are not the seed residues in the path)
     * and maps them to the appropriate position in the topology.
     *
     * @param seedStructure the structure from which to extract match residues
     * @param targetPositions a mapping from residue codes (chain + residue 
     *  number) to indexes in the topology
     * @return a pair consisting of a list of residues from the seed structure,
     *  and the positions in the topology at which these residues should be
     *  added
     */
    pair<vector<Residue *>, vector<int>> getMappedMatchResidues(const Structure &seedStructure, const unordered_map<pair<string, int>, int, pair_hash> &targetPositions);
    int fusePath(const vector<Residue *> &residues, Structure &fusedPath);
    bool pathClashes(const Structure &path, int seedStartIdx);

    /**
     * Finds the next seed to add to the given path segment, and returns the
     * residue adjacent to the segment of the seed that should be added. For
     * example, if appendingForward is true, returns the residue before the
     * part that can be added; if it is false, returns the residue after the
     * part that can be added.
     */
    Residue *findSeedToAdd(const vector<Residue *> &path, int endIndex, const string &currentSeedName, bool appendingForward);
};

#endif /* pathsampler_h */
