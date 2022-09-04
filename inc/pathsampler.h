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
    PathResult(vector<Residue *> originalResidues, Structure fusedPath, vector<vector<Residue*>> topologySeedResidues, int seedStartIdx, fusionOutput fuserScore, int interchainClash = 0, int intrachainClash = 0): _originalResidues(originalResidues), _fusedPath(fusedPath), _topologySeedResidues(topologySeedResidues), _seedStartIdx(seedStartIdx), _fuserScore(fuserScore), _interchainClash(interchainClash), _intrachainClash(intrachainClash) {
        residueSize = _fusedPath.residueSize() - seedStartIdx;
    }
    
    vector<Residue *> getOriginalResidues() { return _originalResidues; }
    int size() { return residueSize; }

    Structure& getFusedStructure() { return _fusedPath; }
    void getFusedPathOnly(Structure &ret) {
        Chain *c = new Chain;
        for (int i = _seedStartIdx; i < _fusedPath.residueSize(); i++) {
            c->appendResidue(new Residue(_fusedPath.getResidue(i)));
        }
        c->setID(_chainID);
        ret.appendChain(c);
    }
    vector<Structure> getTopologySeedResidues() {
        vector<Structure> seedSegments;
        for (int i = 0; i < _topologySeedResidues.size(); i++) {
            seedSegments.emplace_back(_topologySeedResidues[i]);
            Structure& seedSegment = seedSegments[i];
            Residue* pathRes = _originalResidues[i];
            string name = MstUtils::toString(i) + "_" + pathRes->getStructure()->getName() + ":" + MstUtils::toString(pathRes->getNum());
            seedSegment.setName(name);
        }
        return seedSegments;
    }
    fusionOutput getFuserScore() {return _fuserScore;}
    int getIntrachainClash() {return _intrachainClash;}
    int getInterchainClash() {return _interchainClash;}

private:
    vector<Residue *> _originalResidues;
    Structure _fusedPath;
    vector<vector<Residue*>> _topologySeedResidues;
    int _seedStartIdx;
    fusionOutput _fuserScore;
    int _interchainClash;
    int _intrachainClash;
    string _chainID = "0";
    int residueSize = 0; //the residue length of the fused path (no context)
    
};

class PathSampler {
public:
    PathSampler(Structure *target, int overlapLength = 3): overlapLength(overlapLength), _target(target), targetAPV(target->getAtoms()), ps(targetAPV, 10.0) {};
//    ~PathSampler() {};
    
    /**
     Samples the given number of paths from the seed graph or cluster tree.
     */
//    virtual vector<PathResult> sample(int numPaths) = 0;
    
    void setMinimumLength(int length) {minimumLength = length;}
    void setAcceptSingleSeedPaths(bool val) {acceptSingleSeedPaths = val;}
    
    // Optimize without and then with internal coordinate constraints, if set to true
    void setTwoStepFuse(bool val) {twoStepFuse = val;}
    
    static string getPathFromResidues(vector<Residue*> path) {
        string path_string = "";
        for (Residue* R : path) {
            path_string += R->getStructure()->getName() + ":" + MstUtils::toString(R->getResidueIndexInChain());
            if (R != path.back()) path_string += ";";
        }
        return path_string;
    }
    
    /**
     * Generates a PathResult by fusing the given path together, making
     * sure it is valid, and checking for clashes.
     */
    bool emplacePathFromResidues(vector<Residue *> path, vector<PathResult> &results, set<Residue*> fixedResidues = {}, bool ignore_clashes = false);
        
protected:
    Structure *_target = nullptr;
    AtomPointerVector targetAPV;
    ProximitySearch ps;
    
    int overlapLength = 3; // Number of residues to include for overlap when fusing
    int minimumLength = 15;
    bool acceptSingleSeedPaths = false; // If true, accept paths that consist of residues from a single seed
    
    bool twoStepFuse = true;
    
    string chainID = "0";

    AtomPointerVector buildAPV(const vector<Residue *>::const_iterator &begin, const vector<Residue *>::const_iterator &end);

    /**
     * Gets the residues from the target that correspond to the non-seed
     * match residues in the given path.
     *
     * @param pathResidues a list of path residues
     * @return a list of residues from _target that are marked as associated
     *  with the seeds that the path residues come from
     */
    vector<Residue *> getTargetResidues(const vector<Residue *> &pathResidues);
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
    pair<vector<Residue *>, vector<int>> getMappedMatchResidues(const Structure &seedStructure, const map<int,int> &targetPositions);
    int fusePath(const vector<Residue *> &residues, Structure &fusedPath, fusionOutput& fuserScore, vector<vector<Residue *>>& topologySeedResidues, set<Residue*> fixedResidues = {});
    bool pathClashes(const Structure &path, int seedStartIdx, int &interchain_clash, int &intrachain_clash);
};

/**
 Class that randomly samples paths from a set of seed residues which are connected
 by a seed graph.
 */
class SeedGraphPathSampler: public PathSampler {
public:
    /**
     Initialize the path sampler to use a seed graph.
     
     @param target the target structure from which to obtain structure context for fusion
     @param graph a bond-adjacency graph containing seeds to sample from
     */
    SeedGraphPathSampler(Structure *target, SeedGraph *graph, int overlapLength = 3): PathSampler(target, overlapLength), _graph(graph) {
        unordered_set<Residue *> residues = _graph->getResidues();
        _startingResidues = vector<Residue*>(residues.begin(), residues.end());
        if (_startingResidues.size() >= INT_MAX) MstUtils::error("The number of residues in the graph exceeds the max value of an integer.");
    };

    ~SeedGraphPathSampler() {};
    
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
    
    /*
     Prints sampling statistics
     */
    void reportSamplingStatistics() {
        cout << "Total attempts: " << attempts << endl;
        cout << "No overlaps: " << no_overlaps << endl;
        cout << "Redundant: " << redundant << endl;
        cout << "Below minimum length: " << too_short << endl;
        cout << "Contains a clash: " << clashes << endl;
        cout << "Accepted: " << _sampledPaths.size() << endl;
    };
  
    // If true, constrain sampled paths to never use seeds that have been used in previously-sampled paths
    bool uniqueSeeds = false;

    // If non-null, weight possible extensions by the predominant secondary structure classification of the seed
    string *preferredSecondaryStructure = nullptr;
    float secondaryStructureWeight = 2.0f;
  
    // If provided, will always initialize the path from theses residues
    // (only implemented in sampleFromGraph)
    void setStartingPathResidues(vector<Residue*> startingResidues) {
        if (_graph == nullptr) MstUtils::error("Option not available when PathSampler is constructed without a seed graph");
        _initialPathResidues = startingResidues;
        sort(_initialPathResidues.begin(),_initialPathResidues.end(),[](Residue* R1, Residue* R2){return R1->getResidueIndex() < R2->getResidueIndex();});
        cout << "Set the initial path to " << _initialPathResidues.size() << " residues and sorted" << endl;
        if (_initialPathResidues.empty()) MstUtils::error("There are no residues in the provided vector","SeedGraphPathSampler::setStartingResidues");
    }
    
    //Adds the residues of a seed to the fixedResidues set, which is checked before fusing a path
    void addFixedSeed(string seedName) {
        Structure* fixedSeed = _graph->getStructureFromFile(seedName);
        if (fixedSeed == nullptr) MstUtils::error(seedName+" not found in binary file","SeedGraphPathSampler::addFixedSeed");
        vector<Residue*> fixedResidues = fixedSeed->getResidues();
        for (Residue* R : fixedResidues) _fixedResidues.insert(R);
        cout << "Added " << fixedResidues.size() << " fixed residues" << endl;
    }

private:
    SeedGraph *_graph = nullptr;
    
    // The residues that a random residue is sampled from to begin the path
    vector<Residue*> _startingResidues;

    // Used if uniqueSeeds is true
    unordered_set<string> _usedSeeds;
    
    // If not empty, these residues always represent the beginning of a path (sample will
    // extend them in both directions, if possible).
    vector<Residue*> _initialPathResidues;
    
    // If not empty, fusePath will check each path for these residues, and if they are
    // found in the path, those positions will be fixed
    set<Residue*> _fixedResidues;
  
    set<vector<Residue*>> _sampledPaths;
    int attempts = 0, no_overlaps = 0, redundant = 0, too_short= 0, clashes = 0;
    
    vector<Residue*> pathResiduesFromSpecifier(string path_spec);
    
//    int fusePath(const vector<Residue *> &residues, Structure &fusedPath, set<Residue*> fixedResidues, mstreal& fuserScore);
};

/**
 Class that randomly samples paths from a set of seed residues which are connected
 by an overlap cluster tree.
 */
class ClusterTreePathSampler: public PathSampler {
public:
    /**
     Initialize the path sampler to use a cluster tree.
     
     @param target the target structure from which to obtain structure context for fusion
     @param fetcher the fragment fetcher used by the cluster tree
     @param overlapTree an overlap-based cluster tree
     @param overlapSize the number of residues used to determine overlap in the cluster tree
     @param overlapRMSD RMSD cutoff for overlap with a structure to constitute adjacency
     @param residues the complete list of residues from which to sample starting points for paths
     */
    ClusterTreePathSampler(Structure *target, FragmentFetcher *fetcher, ClusterTree *overlapTree, int overlapSize, mstreal overlapRMSD, vector<Residue *> residues): PathSampler(target, overlapSize), _fetcher(fetcher), _searchTree(overlapTree), _searchTreeResidues(residues), overlapRMSD(overlapRMSD) {};

    ~ClusterTreePathSampler() {};
    
    /**
     Samples the given number of paths from the seed graph or cluster tree.
     */
    vector<PathResult> sample(int numPaths);
  
    // If true, constrain sampled paths to never use seeds that have been used in previously-sampled paths
    bool uniqueSeeds = false;

    // If non-null, weight possible extensions by the predominant secondary structure classification of the seed
    string *preferredSecondaryStructure = nullptr;
    float secondaryStructureWeight = 2.0f;
  
private:
    FragmentFetcher *_fetcher = nullptr;
    ClusterTree *_searchTree = nullptr;
    vector<Residue *> _searchTreeResidues;

    // Used if uniqueSeeds is true
    unordered_set<string> _usedSeeds;
  
    // Used only if starting_seed is provided
    vector<Residue*> _startingResidues;
  
    mstreal overlapRMSD = 0.0;

    // Minimum cosine angle between first and last alpha carbons in overlap
    mstreal minCosineAngle = 0.75;

    mstreal cosineAngle(const vector<Residue *> &res1, const vector<Residue *> &res2);

    vector<int> shuffleResultIndexes(ClusterSearchResults &searchResults);

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
