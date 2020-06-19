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

/**
 Object that represents a single path composed of segments from multiple seeds. The
 result contains both 'original residues', which contain parent references to the
 original seeds, a 'fused structure' (both the fused path itself as well as any surrounding
 structural context used to constrain the fusion), and a 'fused path' (only the residues
 in the fused path, with references to original structures removed).
 */
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

/**
 Class that randomly samples paths from a set of seed residues, which are connected
 either by a seed graph or an overlap cluster tree.
 */
class PathSampler {
public:
    /*
     Initialize an empty PathSampler
     */
    PathSampler() {};
    /**
     Initialize the path sampler to use a seed graph.
     
     @param target the target structure from which to obtain structure context for fusion
     @param graph a bond-adjacency graph containing seeds to sample from
     */
    PathSampler(Structure *target, SeedGraph *graph): _target(target), _graph(graph), targetAPV(target->getAtoms()), ps(targetAPV, 10.0) {};

    /**
     Initialize the path sampler to use a cluster tree.
     
     @param target the target structure from which to obtain structure context for fusion
     @param fetcher the fragment fetcher used by the cluster tree
     @param overlapTree an overlap-based cluster tree
     @param overlapSize the number of residues used to determine overlap in the cluster tree
     @param overlapRMSD RMSD cutoff for overlap with a structure to constitute adjacency
     @param residues the complete list of residues from which to sample starting points for paths
     */
    PathSampler(Structure *target, FragmentFetcher *fetcher, ClusterTree *overlapTree, int overlapSize, mstreal overlapRMSD, vector<Residue *> residues): _target(target), _fetcher(fetcher), _searchTree(overlapTree), _searchTreeResidues(residues), overlapLength(overlapSize), overlapRMSD(overlapRMSD), targetAPV(target->getAtoms()), ps(targetAPV, 10.0) {};

    /**
     Samples the given number of paths from the seed graph or cluster tree.
     */
    vector<PathResult> sample(int numPaths);
  
    /**
     Fuse prespecified paths.
     
     @param path_specifiers a list of path specifying strings
     
     A path consisting of residues 1-3 in seed_i and 7-8 in seed_j would be denoted as:
     "seed_i:1,seed_i:2,seed_i:3,seed_j:7,seed_j:8"
     */
  
    vector<PathResult> fusePaths(const vector<string> &path_specifiers);
  
    // If true, constrain sampled paths to never use seeds that have been used in previously-sampled paths
    bool uniqueSeeds = false;

    // If non-null, weight possible extensions by the predominant secondary structure classification of the seed
    string *preferredSecondaryStructure = nullptr;
    float secondaryStructureWeight = 2.0f;
  
    // If provided, will always initialize the path from a residue in this seed
    // (only implemented in sampleFromGraph)
  void setStartingSeed(Structure* seed, string seed_chain) {
      if (_graph == nullptr) MstUtils::error("Option not available when PathSampler is constructed without a seed graph");
    _startingResidues = seed->getChainByID(seed_chain)->getResidues();
    MstUtils::assert((_startingResidues.empty()),"There are no residues in the specified chain");
  }

private:
    Structure *_target = nullptr;
    SeedGraph *_graph = nullptr;
    FragmentFetcher *_fetcher = nullptr;
    ClusterTree *_searchTree = nullptr;
    vector<Residue *> _searchTreeResidues;

    // Used if uniqueSeeds is true
    unordered_set<string> _usedSeeds;
  
    // Used only if starting_seed is provided
    vector<Residue*> _startingResidues;
  
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
  
    vector<Residue*> pathResiduesFromSpecifier(string path_spec);
};

#endif /* pathsampler_h */
