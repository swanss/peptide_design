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

class StructuresBinaryFile {
public:
    StructuresBinaryFile(string filePath, bool read = true): _filePath(filePath), readMode(read) {
        openFileStream(filePath);   
    };

    StructuresBinaryFile(const StructuresBinaryFile &other): readMode(true), _filePath(other._filePath), _filePositions(other._filePositions), _structureNames(other._structureNames) {
        MstUtils::assert(other.readMode, "Copying write-only binary file not supported");
        cout << "Opening file stream for copy, with " << _structureNames.size() << " loaded structure names" << endl;
        openFileStream(_filePath);
    }

    ~StructuresBinaryFile() {
        fs.close();
    }

    Structure *next();
    bool hasNext();
    void skip();
    void reset();
    Structure * getStructureNamed(string name);
    
    void scanFilePositions();

    size_t structureCount() {
        if (_filePositions.empty())
            scanFilePositions();
        return _filePositions.size();
    }

    void jumpToStructureIndex(int idx);

    void insertStructureNames(vector<string> &names) {
        if (_filePositions.empty())
            scanFilePositions();
        names.insert(names.end(), _structureNames.begin(), _structureNames.end());
    }

    void appendStructure(Structure *s);
private:
    string _filePath;
    bool readMode;
    void openFileStream(string filePath);
    fstream fs;
    unordered_map<string, long> _filePositions; // Positions of each structure in the file
    vector<string> _structureNames;
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
     */
    StructureIterator(const vector<string> &filePaths,
                      int batchSize = 1000,
                      vector<string> *chainIDs = nullptr);
    
    StructureIterator(const string binaryFilePath,
                      int batchSize = 1000,
                      string chainID = "0");

    ~StructureIterator() {
        if (binaryFile != nullptr) {
            delete binaryFile;
        }
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

private:
    vector<Structure *> *_lastBatch = nullptr;
    StructuresBinaryFile *binaryFile = nullptr;
    const vector<string> _filePaths;
    vector<string> *_chainIDs;
    string defaultChainID;
    int _batchSize;
    int _batchIndex = 0;
    bool _seedChains;
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
