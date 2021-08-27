#ifndef ClusterUtils_H
#define ClusterUtils_H

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <set>
#include <unordered_map>

#include "msttypes.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "fileutilities.h"
#include "utilities.h"

using namespace MST;

// Abstract parent classes

/**
 * Stores information about a fragment. Subclasses add fields for specific
 * types of fragments.
 */
struct FragmentInfo {
    virtual ~FragmentInfo() = default;
    virtual string toString() = 0;
    virtual FragmentInfo *copy() = 0;
};

/**
 * FragmentFetcher is an abstract base class that handles retrieving fragments
 * from a set of PDB structures. FragmentFetchers store their info in a 
 * FragmentInfo object, which is type-specific.
 */
class FragmentFetcher {
public:
    virtual ~FragmentFetcher() {
        destroyNextResult();   
    };

    /**
     * Checks to see if there is another fragment available to retrieve with
     * the next() method.
     */
    bool hasNext();

    /**
     * Makes and returns the next fragment and fragment info. The client is
     * responsible for memory management of the FragmentInfo pointer. 
     */
    pair<AtomPointerVector, FragmentInfo *> next();
    
    /**
     * Advances to the next fragment without returning it.
     */
    void skip();

    /**
     * Resets the fragment iteration. Subclasses should override this method
     * to define the behavior for their specific iteration technique.
     */
    virtual void reset() = 0;

    /**
     * Retrieves a structure for the given fragment info (which must be an
     * object originally created by this fragment fetcher).
     */
    virtual AtomPointerVector getAPV(FragmentInfo *info) = 0;
    
    /**
     * Defines how many atoms any structure from this fragment fetcher will
     * have.
     */
    virtual int numAtoms() = 0;

    /**
     * Generates the next result to return in next(). Maintains the
     * currentAPVs list for the current structure and advances through it,
     * then fragments the next structure when currentAPVs is exhausted. If
     * there are no more fragments left to yield, sets nextResult to nullptr.
     */
    virtual void makeNextResult() = 0;
    
    /**
     * Destroys the previously-saved value in nextResult, deleting the
     * memory if necessary.
     */
    void destroyNextResult();
    
    /**
     * Initializes a fragment info object from the given string, which should
     * be from the output of a toString call from a FragmentInfo object 
     * of the type associated with this fetcher. The caller is responsible for
     * managing the memory used by the returned object.
     */
    virtual FragmentInfo *makeFragmentInfo(string infoString) = 0;

    virtual Structure *getFullStructure(FragmentInfo *info) = 0;
    virtual Structure getResultStructure(FragmentInfo *info) = 0;
    
    // Helper for caching fragments to return
    pair<AtomPointerVector, FragmentInfo *> *nextResult = nullptr;

    virtual void beginNextSplit() = 0;
};

// FragmentFetcher subclasses 

/**
 * A concrete specialization of FragmentInfo that stores information on a 
 * pair fragment (e.g. two chain fragments sharing a contact). The fragment
 * is defined by a structure index and the atom indexes of the start and end 
 * of each chain, where the atom indexes are numbered from 0 with sidechain
 * atoms removed.
 */
struct PairFragmentInfo: virtual FragmentInfo {
    PairFragmentInfo(int structureIdx, int atomIdx1, int atomIdx2): structureIdx(structureIdx), atomIdx1(atomIdx1), atomIdx2(atomIdx2) {};
    PairFragmentInfo(const PairFragmentInfo &other): structureIdx(other.structureIdx), atomIdx1(other.atomIdx1), atomIdx2(other.atomIdx2) {};
    PairFragmentInfo(string &textRep) {
        vector<string> comps = splitString(textRep, " ");
        MstUtils::assert(comps.size() == 3, "Expected string to have 3 components, got " + to_string(comps.size()));
        structureIdx = stoi(comps[0]);
        atomIdx1 = stoi(comps[1]);
        atomIdx2 = stoi(comps[2]);
    }

    int structureIdx;
    int atomIdx1;
    int atomIdx2;

    string toString() override {
        return to_string(structureIdx) + " " + to_string(atomIdx1) + " " + to_string(atomIdx2);
    }
    FragmentInfo *copy() override {
        return new PairFragmentInfo(structureIdx, atomIdx1, atomIdx2);
    }
};

/**
 * A concrete specialization of FragmentFetcher that works with pair fragments.
 * It uses Fragmenter to break apart each structure using a precompiled list
 * of contacts.
 */
class PairFragmentFetcher: virtual public FragmentFetcher {
public:
    PairFragmentFetcher(string pdbsPath, string contactsPath, int numFlank = 2, int splitCount = 1);
    ~PairFragmentFetcher() override;

    void reset() override;

    AtomPointerVector getAPV(FragmentInfo *info) override;
    int numAtoms() override;

    void makeNextResult() override;

    FragmentInfo *makeFragmentInfo(string infoString) override;
    Structure *getFullStructure(FragmentInfo *info) override;
    Structure getResultStructure(FragmentInfo *info) override;

    void beginNextSplit() override { _splitEnded = false; }
private:
    int _numFlank;
    StructureCache _structureCache;
    vector<string> pdbs;
    vector<string> contactFiles;

    int _splitCount;
    int _splitSize = 0;
    bool _splitEnded = false;

    /**
     * Loads fragments from the structure at the given index into
     * currentAPVs.
     *
     * @return true if the load operation succeeded, false otherwise
     */
    bool loadAPVsFromStructure(int structureIdx);

    vector<AtomPointerVector> currentAPVs;
    int currentStructureIdx = -1;
    Structure *currentStructure = nullptr;
    int currentAPVIdx;
};

struct SingleFragmentInfo: virtual FragmentInfo {
    SingleFragmentInfo(int structureIdx, int atomIdx): structureIdx(structureIdx), atomIdx(atomIdx) {};
    SingleFragmentInfo(const SingleFragmentInfo &other): structureIdx(other.structureIdx), atomIdx(other.atomIdx) {};
    SingleFragmentInfo(string &textRep) {
        vector<string> comps = splitString(textRep, " ");
        MstUtils::assert(comps.size() == 2, "Expected string to have 2 components, got " + to_string(comps.size()));
        structureIdx = stoi(comps[0]);
        atomIdx = stoi(comps[1]);
    }

    int structureIdx;
    int atomIdx;

    string toString() override {
        return to_string(structureIdx) + " " + to_string(atomIdx);
    }
    FragmentInfo *copy() override {
        return new SingleFragmentInfo(structureIdx, atomIdx);
    }

};

class SingleFragmentFetcher: virtual public FragmentFetcher {
public:
    SingleFragmentFetcher(StructuresBinaryFile *binaryFile, int numResidues, string chainID = "", int splitCount = 1);
    SingleFragmentFetcher(StructuresBinaryFile *binaryFile, StructureCache *cache, int numResidues, string chainID = "", int splitCount = 1);
    //SingleFragmentFetcher(string binaryFilePath, int numResidues, string chainID = "");
    ~SingleFragmentFetcher();
    void reset() override;

    AtomPointerVector getAPV(FragmentInfo *info) override;
    int numAtoms() override;

    void makeNextResult() override;
    FragmentInfo *makeFragmentInfo(string infoString) override;
    Structure *getFullStructure(FragmentInfo *info) override;
    Structure getResultStructure(FragmentInfo *info) override;

    /// If set, only yield residues from the chain with this ID
    string chainID = "";

    vector<Residue *> getAllResidues() {
        vector<string> fileNames;
        _binaryFile->insertStructureNames(fileNames);
        vector<Residue *> residues;
        for (string name: fileNames) {
            Structure *s = _structureCache->getStructure(name);
            vector<Residue *> structureRes;
            if (chainID.empty())
                structureRes = s->getResidues();
            else {
                Chain *c = s->getChainByID(chainID);
                if (c)
                    structureRes = c->getResidues();
            }
            residues.insert(residues.end(), structureRes.begin(), structureRes.end());
        }
        return residues;
    }

    /// Allows the fragment fetcher to start returning fragments again
    void beginNextSplit() override { _splitEnded = false; }
    
private:
    int _numResidues;
    vector<string> pdbs;
    StructuresBinaryFile *_binaryFile = nullptr;
    StructuresBinaryFile *_cacheBinaryFile = nullptr;
    StructureCache *_structureCache = nullptr;

    bool ownsStructureCache = true;

    int _splitCount;
    int _splitSize = 0;
    bool _splitEnded = false;

    /**
     * Loads fragments from the structure at the given index into
     * currentAPVs.
     *
     * @return true if the load operation succeeded, false otherwise
     */
    bool loadAPVsFromStructure(int structureIdx);

    vector<AtomPointerVector> currentAPVs;
    int currentStructureIdx = -1;
    Structure *currentStructure = nullptr;
    int currentAPVIdx = -1;
};

/**
 * A specialization of FragmentFetcher that wraps another FragmentFetcher,
 * but returns only one in every N fragments (where N is the number of 
 * workers). The worker index should be zero-indexed.
 */
class BatchWorkerFragmentFetcher: virtual public FragmentFetcher {
public:
    BatchWorkerFragmentFetcher(FragmentFetcher *fetcher, int numWorkers, int workerIndex): _fetcher(fetcher), _numWorkers(numWorkers), _workerIndex(workerIndex) {};
    void reset() override;

    AtomPointerVector getAPV(FragmentInfo *info) override;
    int numAtoms() override;

    void makeNextResult() override;
    FragmentInfo *makeFragmentInfo(string infoString) override;
    Structure *getFullStructure(FragmentInfo *info) override;
    Structure getResultStructure(FragmentInfo *info) override;

    void beginNextSplit() override { _fetcher->beginNextSplit(); }
private:
    FragmentFetcher *_fetcher;
    int _numWorkers;
    int _workerIndex;
    int _fragmentIndex = 0;
};

/**
 * A FragmentFetcher that enables the live creation of fragments. It does
 * not support iteration, but rather allows you to add fragments one at a
 * time, then call insert() on the cluster tree to add them yourself.
 *
 * OnlineFragmentFetcher copies the FragmentInfos that you pass to it.
 */
class OnlineFragmentFetcher: virtual public FragmentFetcher {
public:
    OnlineFragmentFetcher() {};
    void reset() override {};

    AtomPointerVector getAPV(FragmentInfo *info) override { return fragments[info]; };
    int numAtoms() override { return _numAtoms; };

    void makeNextResult() override {};
    FragmentInfo *makeFragmentInfo(string infoString) override {
        MstUtils::error("OnlineFragmentFetcher creation of fragment infos not supported");   
        return nullptr;
    };
    Structure *getFullStructure(FragmentInfo *info) override {
        MstUtils::error("OnlineFragmentFetcher full structure retrieval not supported");   
        return nullptr;
    };
    Structure getResultStructure(FragmentInfo *info) override {
        MstUtils::error("OnlineFragmentFetcher structure retrieval not supported");   
        return Structure();
    };

    void insert(const AtomPointerVector &apv, FragmentInfo *info) {
        if (_numAtoms == 0) {
            _numAtoms = apv.size();
        } else {
            MstUtils::assert(apv.size() == _numAtoms, "APV must be " + to_string(_numAtoms) + " atoms, not " + to_string(apv.size()));
        }
        fragments[info->copy()] = apv;
    }

private:
    unordered_map<FragmentInfo *, AtomPointerVector> fragments;
    int _numAtoms = 0;
};

// Other helpers


/**
 * A helper class that stores previously computed RMSDs.
 */
class RMSDCache {
public:
    RMSDCache(FragmentFetcher *fetcher, bool sharedCoordinates = false): _fetcher(fetcher), _sharedCoordinates(sharedCoordinates) {};

    double bestRMSD(FragmentInfo *f1, FragmentInfo *f2) {
        if (_cache.count(f1) != 0 && _cache[f1].count(f2) != 0)
            return _cache[f1][f2];
        else if (_cache.count(f2) != 0 && _cache[f2].count(f1) != 0)        
            return _cache[f2][f1];
        else if (f1 == f2)
            return 0.0;
        AtomPointerVector apv1 = _fetcher->getAPV(f1);
        AtomPointerVector apv2 = _fetcher->getAPV(f2);
        double ret = _sharedCoordinates ? _calc.rmsd(apv1, apv2) : _calc.bestRMSD(apv1, apv2);
        _cache[f1][f2] = ret;
        return ret;
    }

    double bestRMSD(FragmentInfo *f1, const vector<Atom *> &apv2) {
        return _sharedCoordinates ? _calc.rmsd(_fetcher->getAPV(f1), apv2) : _calc.bestRMSD(_fetcher->getAPV(f1), apv2);
    }

    double bestRMSD(const vector<Atom *> &apv1, const vector<Atom *> &apv2) {
        return _sharedCoordinates ? _calc.rmsd(apv1, apv2) : _calc.bestRMSD(apv1, apv2);
    }

    void setQuery(AtomPointerVector const* query) {
        _queryCache.clear();
        _query = query;
    }

    double bestQueryRMSD(FragmentInfo *f) {
        if (_queryCache.count(f) == 0)
            _queryCache[f] = _sharedCoordinates ? _calc.rmsd(_fetcher->getAPV(f), *_query) : _calc.bestRMSD(_fetcher->getAPV(f), *_query);
        return _queryCache[f];
    }
private:
    FragmentFetcher *_fetcher = nullptr;
    RMSDCalculator _calc;
    unordered_map<FragmentInfo *, unordered_map<FragmentInfo *, double>> _cache;
    unordered_map<FragmentInfo *, double> _queryCache;
    AtomPointerVector const* _query = nullptr;
    bool _sharedCoordinates = false;
};

// Greedy Clustering Utils
struct seedWindowInfo {
public:
    seedWindowInfo() {};
    seedWindowInfo(string _structure_name, int _position) : structure_name(_structure_name), position(_position) {};
    
    string getName() {
        return structure_name + "_" + MstUtils::toString(position);
    }
    
    bool operator<(const seedWindowInfo& other) const {
        if (structure_name != other.structure_name) {
            return structure_name < other.structure_name;
        } else {
            return position < other.position;
        }
    }
    
    string structure_name;
    int position;
};

class GreedyClusterer {
public:
    GreedyClusterer(string seedBin_path, int window_size);
    
    ~GreedyClusterer() {
        delete seedBin;
    }
    
    void addOverlapInfo(FuseCandidateFile file);
    
    // procedure borrowed from MST/src/msttypes.cpp Clusterer::greedyCluster
    void performClustering(mstreal max_coverage);
    
    void writeClusterInfo(string path, Chain* peptide = nullptr, bool pdbs = true, bool verbose = false);
    
    Structure getSeedWindowStructure(seedWindowInfo info) {
        Structure* extendedfragment = seedStructures.getStructure(info.structure_name);
        vector<Residue*> seed_residues = extendedfragment->getChainByID(seed_chain_id)->getResidues();
        vector<Residue*> seed_window_residues(seed_residues.begin()+info.position,seed_residues.begin()+info.position+window_length);
        return (seed_window_residues);
    }
    
private:
    int window_length;
//    mstreal cluster_radius = 1.0;
//    int max_cluster_number = 10000;
//    mstreal max_coverage = 1.0;
    string seed_chain_id = "0";
    
    /*
     The key is a unique seed window ID and the value is the information needed to retrieve the seed
     window from the structure cache.
     */
    map<int,seedWindowInfo> seedWindows;
    map<seedWindowInfo,int> seedWindowsRev; //to lookup ID by seed name + position
    
    /*
     The key is the seed window ID and the mapped value is all overlapping seed windows. For convenience,
     each overlap is listed twice in the map. This means that if seed window A and B have an overlap, this can
     be accessed by querying the map with the ID of either A or B.
     */
    map<int,set<int>> seedWindowOverlaps;
    
    /*
     Keeps track of the seed windows that have not yet been covered in a cluster
     */
    set<int> remainingSeedWindows;
    int numTotalSeedWindows;
    
    /* Stores all of the cluster data
     pair<cluster_ID,set<cluster_ID>>
     */
    vector<pair<int,set<int>>> clusters;
    
    StructuresBinaryFile* seedBin;
    StructureCache seedStructures;
    
};

#endif
