//
//  fileformats.h
//  dtermen
//
//  Created by Venkatesh Sivaraman on 1/28/19.
//

#ifndef fileformats_h
#define fileformats_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <set>
#include <map>
#include <list>
#include <unordered_map>
#include <fstream>

#include "msttypes.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "structure_iter.h"

using namespace std;
using namespace MST;

vector<string> splitString(string s, const string &delim);

/**
 Determines the target residue from which the given seed came. Assumes the file
 name is of the form "XXXX_XABC_match_XX_seed_X.pdb", and extracts the number
 "ABC".
 
 @param seedName the file name of the seed
 @return the target residue index for this seed, or -1 if one was not found
 */
int getTargetResidueIndex(string seedName);

/**
 * Determines the target residue code, including the chain ID and the number,
 * assuming the file name is of the form "[TARGETNAME]-[CHAINID][NUMBER]-etc.pdb".
 */
pair<string, int> getTargetResidueCode(string seedName);

/**
 Type that encapsulates a CSV file indicating paths for seeds. At minimum, the
 file should contain the path for each seed and a string indicating which chains
 can be used as seeds (e.g. "AB" for both chains A and B). It may also contain
 one or more scores for these seeds.
 
 The CSV file should be formatted with a header line, where the first column is
 called "path" and the second column is called "chids".
 */
class SeedListFile {
public:
    SeedListFile(string path): _path(path) {}
    ~SeedListFile();
    
    /**
     Reads the file.
     
     @param pdbPrefix a path to prepend to each path in the file for return
     @return a list of seed paths from the file, and a list of the chain IDs
     corresponding to each seed path
     */
    pair<vector<string>, vector<string>> read(string pdbPrefix = "");
    /**
     @return the names of each metadata field in the file
     */
    vector<string> getMetadataFields() { return _metadataNames; }
    /**
     @return the stored metadata for the given path
     */
    vector<string> getMetadata(string path) { return metadata[path]; }
    /**
     Sets the metadata field names for writing in the CSV file.
     
     @param metadataNames the new field names
     */
    void setMetadataFields(vector<string> metadataNames) { _metadataNames = metadataNames; }
    /**
     Writes the given seed path to the file, with the optional scores provided.
     
     @param seed the seed path to write
     @param chainIDs a concatenated string of the seed chain IDs
     @param metadata a list of metadata items to provide in the file (must have
     previously provided a set of headers of the same length)
     */
    void write(string seed, string chainIDs, vector<string> metadata = vector<string>());
private:
    vector<string> _metadataNames;
    vector<string> listedPaths;
    vector<string> listedChainIDs;
    map<string, vector<string>> metadata; // keys paths to metadata fields
    string _path;
    ifstream *readstream = nullptr;
    ofstream *writestream = nullptr;
};

struct FuseCandidate {
    string file1;
    string chain1;
    string file2;
    string chain2;
    Structure *structure1 = nullptr;
    Structure *structure2 = nullptr;
    
    int overlapPosition1;
    int overlapPosition2;
    int overlapSize;
    float rmsd;
    
    /**
     Sets the properties of the FuseCandidate based on the given structure.
     @param structure the structure to use as the first structure in the candidate
     @param chainID the ID of the seed chain within the structure
     */
    void setStructure1(Structure *structure, string chainID) {
        structure1 = structure;
        file1 = structure->getName();
        chain1 = chainID;
    }
    /**
     Sets the properties of the FuseCandidate based on the given structure.
     @param structure the structure to use as the second structure in the candidate
     @param chainID the ID of the seed chain within the structure
     */
    void setStructure2(Structure *structure, string chainID) {
        structure2 = structure;
        file2 = structure->getName();
        chain2 = chainID;
    }
    
    /**
     Loads the structure1 and structure2 properties of the fuse candidate from
     the candidate's metadata. The user is responsible for freeing the structures
     using freeStructures().
     @param pathPrefix a path to prepend to the file1 and file2 attributes to
     find the PDB path
     */
    void loadStructures(string pathPrefix = "") {
        structure1 = new Structure(pathPrefix + file1);
        structure2 = new Structure(pathPrefix + file2);
    }
    
    /**
     Deletes the structures associated with this fuse candidate.
     */
    void freeStructures() {
        delete structure1;
        delete structure2;
    }
};

/**
 Type that encapsulates a CSV file listing fuse candidates, allowing for read
 and write operations. Note: A FuseCandidateFile can be used for reading or
 writing, but not both.
 */
class FuseCandidateFile {
public:
    FuseCandidateFile(string path, bool shouldAppend = false): _path(path), _shouldAppend(shouldAppend) {}
    ~FuseCandidateFile();
    
    /**
     Reads up to numLines lines from the file. May return fewer if the end of
     the file is reached.
     
     @param numLines the number of lines to read
     @return a list of fuse candidates from the file, with up to numLines entries
     */
    vector<FuseCandidate> read(int numLines);
    /**
     Writes the given list of fuse candidates to the file.
     
     @param candidates the list of candidates to write
     @param parentPath the path to write other paths relative to
     */
    void write(vector<FuseCandidate> candidates, string parentPath = "");
    /**
     Writes the given fuse candidate to the file.
     
     @param candidate the candidate to write
     @param parentPath the path to write other paths relative to
     */
    void write(FuseCandidate candidate, string parentPath = "");
private:
    bool _isFileDone;
    string _path;
    bool _shouldAppend;
    ifstream *readstream = nullptr;
    ofstream *writestream = nullptr;
};

#define MAX_COORDINATE 1000

/**
 Represents a series of labeled points in 3D space. All points are assumed to
 lie within MAX_COORDINATE of the origin along all axes. Note that only the
 labels and their rough locations are stored in the map.
 */
template <class T>
class PositionMap {
public:
    /**
     Initializes a position map where the hash bins are separated by interval.
     
     @param interval the side length of each cubic hash bin
     */
    PositionMap(int interval = 5): _interval(interval), _maxCoordBin(MAX_COORDINATE / interval * 2) {
        _maxBin = _maxCoordBin * _maxCoordBin * _maxCoordBin;
    }
    
    /**
     Insert the labeled point in the bin for the given location.
     
     @param x the x position
     @param y the y position
     @param z the z position
     @param identifier the label to insert
     */
    void insert(float x, float y, float z, T identifier);
    /**
     Insert the labeled point in every bin that intersects with the rectangular
     prism bounded by the given min and max coordinates.
     
     @param xmin the minimum x position
     @param ymin the minimum y position
     @param zmin the minimum z position
     @param xmax the maximum x position
     @param ymax the maximum y position
     @param zmax the maximum z position
     @param identifier the label to insert
     */
    void insert(float xmin, float ymin, float zmin, float xmax, float ymax, float zmax, T identifier);
    /**
     Finds labels near a given position.
     
     @param x the x position near which to search
     @param y the y position near which to search
     @param z the z position near which to search
     @return a list of labels whose locations are near the given position
     */
    vector<T> itemsNear(float x, float y, float z);
private:
    int _interval;
    long _maxCoordBin, _maxBin;
    unordered_map<long, vector<T>> _table;
    
    long hash(float x, float y, float z);
    vector<long> hash(float xmin, float ymin, float zmin, float xmax, float ymax, float zmax);
    vector<long> nearbyBins(float x, float y, float z);
    void pushIfInRange(long bin, vector<long> &arr);
};

/**
 StructureCache maintains a shared set of Structure objects, which provides
 convenient lookup and iteration. This is useful because keeping the same set of
 Structures allows one to lookup residues by their pointers (otherwise, Residue
 objects would be duplicated but represent the same residue).
 
Note: StructureCache has a capacity parameter that saves memory by using a least-
recently used (LRU) caching algorithm to save frequently-used structures. However,
this means that structure data may get deleted from the cache, including any
auxiliary objects owned by the structure. To avoid segmentation faults, copy
any structure data that needs to be persisted over many calls to the structure
cache, or omit the capacity argument to use an effectively infinite-sized cache.
 */
class StructureCache {
public:
    /**
     Initializes a structure cache that will load files from the given directory,
     i.e. the pdb prefix will be prepended to any paths that are requested.
     */
    StructureCache(string pdbPrefix = "", long capacity = 100000000): pdbPrefix(pdbPrefix), capacity(capacity) { };
    
    StructureCache(StructuresBinaryFile *binaryFile, long capacity = 100000000): binaryFile(binaryFile), capacity(capacity) { };

    ~StructureCache();
    
    StructureCache(const StructureCache& other): pdbPrefix(other.pdbPrefix), cache(other.cache), binaryFile(other.binaryFile), cachePointers(other.cachePointers), capacity(other.capacity) { };
    
    /// If true, cache a vector of pointers to each structure's atoms.
    bool storeAtoms = false;

    /**
     * Preloads all structures from the binary file, up to
     * the capacity of the cache.
     */
    void preloadFromBinaryFile();

    /**
     Returns the structure for the given name/path, loading it fresh if it is not
     already loaded. Uses the given path prefix if provided, otherwise uses the
     cache's default one (provided in the constructor).
     
     @param name the name of the PDB file
     @param prefix a path to prepend the file name with (if empty, uses the
            StructureCache's default)
     @return a Structure loaded from the PDB file
     */
    Structure *getStructure(string name, string prefix = "");
    
    /**
     * Performs the same action as getStructure, but returns a vector of
     * atoms. The storeAtoms property must be true for this to work.
     */
    vector<Atom *> &getAtoms(string name, string prefix = "");

    /**
     Determine whether this cache already has the given structure loaded.
     
     @param name the name of the PDB file
     @return true if the structure is already loaded, and false if not
     */
    bool hasStructure(string name);
    
    /**
     Removes the structure from the cache if it is present.

     @param name the name of the PDB file
    */
    void removeStructure(string name);

    /**
     * Removes all structures from the cache.
     */
    void clear();

    /** @return the PDB path prefix used by this cache */
    string getPDBPrefix() { return pdbPrefix; }
    
    /** @param prefix the new prefix to use */
    void setPDBPrefix(string prefix) { pdbPrefix = prefix; }
    
    /**
     Allows clients to iterate over the structures in this cache (in no particular
     order).
     */
    class iterator {
    public:
        typedef iterator self_type;
        typedef Structure * value_type;
        typedef int difference_type;
        typedef forward_iterator_tag iterator_category;
        iterator(list<Structure *>::iterator it) : it_(it) { }
        iterator(const iterator& other): it_(other.it_) { }
        self_type operator++() { it_++; return *this; }
        self_type operator++(int junk) { self_type i = *this; it_++; return i; }
        Structure * operator*() { return *it_; }
        Structure * const * operator->() { return &(*it_); }
        bool operator==(const self_type& rhs) { return it_ == rhs.it_; }
        bool operator!=(const self_type& rhs) { return it_ != rhs.it_; }
    private:
        list<Structure *>::iterator it_;
    };
    
    /**
     Iterator pointing to the beginning of the structure cache.
     */
    iterator begin() {
        iterator it(cache.begin());
        return it;
    }
    
    /**
     Iterator pointing to the end of the structure cache.
     */
    iterator end() {
        iterator it(cache.end());
        return it;
    }
    
private:
    string pdbPrefix;
    long capacity;
    unordered_map<string, list<Structure *>::iterator> cachePointers;
    unordered_map<string, vector<Atom *>> atoms;
    list<Structure *> cache;
    StructuresBinaryFile *binaryFile = nullptr;
};

#endif /* fileformats_h */
