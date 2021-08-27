//
//  structure_iter.h
//  dtermen
//
//  Created by Venkatesh Sivaraman on 12/19/18.
//

#ifndef structure_iter_h
#define structure_iter_h

#include <stdio.h>
#include "msttypes.h"
#include <unordered_map>
#include <list>

using namespace std;
using namespace MST;

/**
 Utility function that extracts a substructure consisting of the given chains.
 
 @param s the structure from which to extract
 @param chainIDs the chain IDs to extract
 @param newS the destination structure for the new chains
 */
void extractChains(Structure &s, vector<string> chainIDs, Structure &newS);
void extractChains(Structure &s, string chainIDs, Structure &newS);


/*
 StructuresBinaryFile handles reading and writing of binary files containing seeds.
 
 In read mode, the class automatically detects version and allows the user to read structures
 one-by-one. If the file is version 2, it can also contain meta-data pertaining to each of the seeds.
 This can be accessed by first calling scanFilePositions() to read all meta-data in the file, and then
 getStructureProperty() to get specific values.
 
 In write mode, only version 2 is supported. By default, a new file will be opened at the specified
 path and any existing file will be overwritten. However, if append is specified, the new structures
 will be appended to the existing file
 
 Note: seed binary files written before 06/28/20 are all version 1, and all files written afterwards are
 version 2.
 */

class StructuresBinaryFile {
public:
    StructuresBinaryFile(string filePath, bool read = true, bool append = false): _filePath(filePath), readMode(read), _append(append){
        detectFileVersion();
        cout << "read mode: " << readMode << "\tappend: " << _append << "\tversion: " << _version << "\topening file: " << filePath << endl;
        openFileStream();
    };

    StructuresBinaryFile(const StructuresBinaryFile &other): readMode(true), _filePath(other._filePath), _filePositions(other._filePositions), _structureNames(other._structureNames), _version(other._version) {
        MstUtils::assert(other.readMode, "Copying write-only binary file not supported");
        cout << "Opening file stream for copy, with " << _structureNames.size() << " loaded structure names" << endl;
        openFileStream();
    }

    ~StructuresBinaryFile() {
        if ((!readMode) && (structure_added)) MstUtils::writeBin(fs,'E');
        fs.close();
    }

    Structure *next();
    bool hasNext();
    void skip();
    void reset();
    Structure * getStructureNamed(string name);
    
    /*
     Finds the position of each structure section in the file and populates the structure
     property variables.
     */
    void scanFilePositions();

    size_t structureCount() {
        if (_filePositions.empty())
            scanFilePositions();
            reset();
        return _filePositions.size();
    }
    
    int getStructurePropertyInt(string prop, string name) {return seed_dscrt_vals.at(prop).at(name);}
    mstreal getStructurePropertyReal(string prop, string name) {return seed_real_vals.at(prop).at(name);}

    void jumpToStructureIndex(int idx);

    void insertStructureNames(vector<string> &names) {
        if (_filePositions.empty())
            scanFilePositions();
        names.insert(names.end(), _structureNames.begin(), _structureNames.end());
    }
    
    //Only version 1 is supported for writing
    void appendStructure(Structure *s);
    
    //All properties are appended to the section following the last added structure
    void appendStructurePropertyInt(string prop, int val);
    void appendStructurePropertyReal(string prop, mstreal val);
  
    set<string> getPropertyNamesInt() {
        set<string> properties;
        for (auto element : seed_dscrt_vals) properties.insert(element.first);
        return properties;
      
    }
    set<string> getPropertyNamesReal() {
        set<string> properties;
        for (auto element : seed_real_vals) properties.insert(element.first);
        return properties;
    }
    
    vector<string> getStructureNames();
    
private:
    string _filePath;
    bool readMode;
    /*
     0: the original seed binary format
     1: updated version with meta-data (06/21/20)
     */
    int _version;
    bool _append;
    bool structure_added = false;
    
    void openFileStream();
    void detectFileVersion();
    pair<Structure*,long> readNextFileSection(bool save_metadata = false);

    fstream fs;
    long fileStart;
    unordered_map<string, long> _filePositions; // Positions of each structure in the file
    vector<string> _structureNames;
    
    /*
     The following variables are only populated if scanFilePositions is called for a version 1
     seed binary file.
     
     seed_vals[property_name][structure_name] = value;
     */
    map<string,map<string,int>> seed_dscrt_vals;
    map<string,map<string,mstreal>> seed_real_vals;
};

/**
 StructureCache maintains a shared set of Structure objects, which provides
 convenient lookup and iteration. This is useful because keeping the same set of
 Structures allows one to lookup residues by their pointers (otherwise, Residue
 objects would be duplicated but represent the same residue).
 
Note: StructureCache has a capacity parameter that saves memory by using a least-
recently used (LRU) caching algorithm to save frequently-used structures. However,
this means that if the capacity is exceeded, structure data will be deleted from
the cache, including any auxiliary objects owned by the structure. To avoid
segmentation faults, copy any structure data that needs to be persisted over many
calls, or omit the capacity argument to use an effectively infinite-sized cache.
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
     * Preloads all structures from binary file, up to
     * the capacity of the cache.
     */
    void preloadFromBinaryFile();
    
    void preloadFromPDBList(string pdbList);
    
    bool isPreloaded() {return preloaded;}
    bool belowCapacity() {return _belowCapacity;}
    int size() {return preloaded ? cache.size() : -1;}
    
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
//        //this method results in the iterator skipping end()
//        self_type advance(int step_size) {
//            self_type i = *this;
//            for (int i = 0; i < step_size; i++) {
//                it_++;
//                if (it_ == it_.end()) break;
//            }
//            return i;
//        }
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
    bool preloaded = false;
    bool _belowCapacity = false;
    unordered_map<string, list<Structure *>::iterator> cachePointers;
    unordered_map<string, vector<Atom *>> atoms;
    list<Structure *> cache;
    StructuresBinaryFile *binaryFile = nullptr;
};


class PairStructureIterator;

/**
 Helper class that manages a list of PDB structure paths and loads their
 contents in batches.
 */
class StructureIterator {
    friend class PairStructureIterator;
public:
    /**
     @param filePaths the file paths to load
     @param batchSize the number of structures to load in each iteration
     @param chainIDs the chain IDs to use from the file paths; if non-null,
            should be the same length as filePaths
     @param workerIndex assigns the worker to specific batches, should be in [1,numWorkers]
     @param numWorkers the total number of workers
     */
    StructureIterator(const vector<string> &filePaths,
                      int batchSize = 1000,
                      vector<string> *chainIDs = nullptr,
                      int workerIndex = 1,
                      int numWorkers = 1);
    
    StructureIterator(const string binaryFilePath,
                      int batchSize = 1000,
                      string chainID = "0",
                      int workerIndex = 1,
                      int numWorkers = 1);

    ~StructureIterator() {
        if (binaryFile != nullptr) {
            delete binaryFile;
        }
        if (_currentBatch != nullptr)
          delete _currentBatch;
    }
    /**
     @return whether there are structures that have not been loaded and returned
             yet
     */
    bool hasNext();
    /**
     Resets the iterator to start at the beginning of its file path list.
     */
    void reset();
    /**
     Loads and returns the next batch of structures. hasNext() must be true
     before calling this method.
     
     @return the next batch of structures
     */
    vector<Structure *> next();
    /**
     Skips the current batch of structures.
     */
    void skip();
  
    int getBatchSize() {return _batchSize;}
    
    /**
     Sets whether the class has ownership of the seed structures (is responsible for deleting them) or not.
     */
    void hasOwnership(bool val) {maintainsOwnership = val;}

private:
    vector<Structure *> *_currentBatch = nullptr;
    StructuresBinaryFile *binaryFile = nullptr;
    const vector<string> _filePaths;
    vector<string> *_chainIDs;
    string defaultChainID;
    int _batchSize;
    int _batchIndex;
    bool _seedChains;
    int _workerIndex = 0;
    int _numWorkers = 1;
    bool maintainsOwnership = true;
    
    bool nextResultAvailable = false;
    
    void makeNextResult();
};

/**
 Helper class that manages a list of PDB structure files and returns their
 batches in pairwise combinations. For example, to iterate over pairs of
 structures from the list
 
 A, B, C, D, E, F, G, H, I,
 
 in batches of size 3, the resulting batch pairs would be (in some order)
 
 <A, B, C>, <A, B, C>
 <A, B, C>, <D, E, F>
 <A, B, C>, <G, H, I>
 <D, E, F>, <A, B, C>
 <D, E, F>, <D, E, F>
 <D, E, F>, <G, H, I>
 <G, H, I>, <A, B, C>
 <G, H, I>, <D, E, F>
 <G, H, I>, <G, H, I>
 
 As can be seen from this output, you would still have to iterate over pairwise
 combinations of the two lists in order to test the full Cartesian product.
 */
class PairStructureIterator {
public:
    /**
     @param filePaths the file paths to load
     @param batchSize the number of structures to load in each iteration
     @param chainIDs the chain IDs to use from the file paths; if non-null,
            should be the same length as filePaths
     */
    PairStructureIterator(const vector<string> &filePaths,
                          int batchSize = 1000,
                          vector<string> *chainIDs = nullptr);
    
    PairStructureIterator(const string binaryFilePath,
                          int batchSize = 1000,
                          string chainID = "0");

    /**
     * Indicates whether to remove half the batches that are simply
     * the reverse order of other batches.
     */
    bool skipSymmetricBatches = true;
    /**
     @return whether there are structures that have not been loaded and returned
     yet
     */
    bool hasNext();
    /**
     Resets the iterator to start at the beginning of its file path list.
     */
    void reset();
    /**
     Loads and returns the next batch pair of structures. hasNext() must be true
     before calling this method.
     
     @return the next batch pair of structures
     */
    tuple<vector<Structure *>, vector<Structure *>> next();
    /**
     Skips the current batch of structures.
     */
    void skip();
private:
    vector<Structure *> _currentFirst;
    StructureIterator _first;
    StructureIterator _second;
};

class BatchPairStructureIterator {
public:
    BatchPairStructureIterator(const string &binaryFilePath,
                               int workerIndex = 0,
                               int numWorkers = 1,
                               int batchSize = 1000,
                               string chainID = "0");

    ~BatchPairStructureIterator();

    bool hasNext();
    void reset();
    pair<vector<Structure *>, vector<Structure *>> next();
    int getFirstIndex() { return firstIndex; }
    int getSecondIndex() { return secondIndex; }
    
private:
    StructuresBinaryFile *binaryFile = nullptr;
    int workerIndex;
    int numWorkers;
    int batchSize;
    int firstIndex = -1;
    int secondIndex = -1;
    int numRows = 0; // number of batches in a single loop through the structures
    string chainID;

    vector<Structure *> currentFirst;
    vector<Structure *> currentSecond;
    bool nextResultAvailable = false;

    void makeNextResult();
};
#endif /* structure_iter_h */
