//
//  structure_iter.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 12/19/18.
//

#include <stdio.h>
#include "structure_iter.h"
#include "mstsystem_exts.h"
#include <regex>
#include <vector>

using namespace std;
using namespace MST;

vector<string> splitChainIDs(const string &fn) {
    vector<string> result;
    for (int i = 0; i < fn.length(); i++) {
        result.push_back(fn.substr(i, 1));
    }
    
    return result;
}

void extractChains(Structure &s, vector<string> chainIDs, Structure &newS) {
    //    cout << "chains: ";
    //    for (auto& val: chainIDs) cout << val << ", ";
    //    cout << endl;
    
    if (chainIDs.size() == 0) {
        cerr << "no chains to create substructure!" << endl;
        newS = s;
        return;
    }
    newS.setName(s.getName());
    for (auto cID: chainIDs) {
        Chain *theChain = new Chain(*s.getChainByID(cID));
        newS.appendChain(theChain);
    }
}

void extractChains(Structure &s, string chainIDs, Structure &newS) {
    extractChains(s, splitChainIDs(chainIDs), newS);
}

vector<string> seedChainIDs(const string &fn) {
    regex chainIDRegex("chids([A-Za-z0-9]+)_");
    std::smatch base_match;
    vector<string> result;
    if (!std::regex_search(fn, base_match, chainIDRegex)) {
        cerr << fn << endl;
        return result;
    }
    
    if (base_match.size() != 2) {
        cerr << "match wrong size: " << fn << endl;
        return result;
    }
    std::ssub_match base_sub_match = base_match[1];
    std::string base = base_sub_match.str();
    for (int i = 0; i < base.length(); i++) {
        result.push_back(base.substr(i, 1));
    }
    
    if (result.size() == 0)
        cerr << "Empty: " << fn << endl;
    return result;
}

void StructuresBinaryFile::openFileStream(string filePath) {
    if (readMode)
        MstUtils::openFile(fs, filePath, fstream::in | fstream::binary, "StructuresBinaryFile::openFileStream");
    else
        MstUtils::openFile(fs, filePath, fstream::out | fstream::binary | fstream::app, "StructuresBinaryFile::openFileStream");
}

bool StructuresBinaryFile::hasNext() {
    MstUtils::assert(readMode, "hasNext not supported in write mode");
    return fs.peek() != EOF;
}

Structure *StructuresBinaryFile::next() {
    MstUtils::assert(readMode, "next not supported in write mode");
    Structure* S = new Structure();
    S->readData(fs);
    return S;
}

void StructuresBinaryFile::skip() {
    MstUtils::assert(readMode, "skip not supported in write mode");
    Structure *s = new Structure();
    s->readData(fs);
    delete s;
}

void StructuresBinaryFile::reset() {
    MstUtils::assert(readMode, "reset not supported in write mode");
    cout << "Resetting" << endl;
    fs.clear(); // this is necessary in case ifstream doesn't clear eofbit
    fs.seekg(0, fs.beg);
}

void StructuresBinaryFile::scanFilePositions() {
    //MstUtils::assert(readMode, "scanFilePositions not supported in write mode");
    if (!_structureNames.empty())
        return;
    cout << "Scanning file positions..." << endl;
    fs.seekg(0, fs.beg);
    _structureNames.clear();
    while (fs.peek() != EOF) {
        Structure *S = new Structure();
        long pos = fs.tellg();
        S->readData(fs);
        if (pos < 0) {
            cout << pos << endl;
        }
        _filePositions[S->getName()] = pos;
        _structureNames.push_back(S->getName());
        delete S;
    }
    cout << "Done scanning file" << endl;
}

Structure * StructuresBinaryFile::getStructureNamed(string name) {
    MstUtils::assert(readMode, "getStructureNamed not supported in write mode");
    if (_filePositions.size() == 0) {
        scanFilePositions();
    }
    if (_filePositions.count(name) == 0) {
        cout << "Structure doesn't exist" << endl;
        return nullptr;
    }
    fs.seekg(_filePositions[name], fs.beg);
    Structure *S = new Structure();
    S->readData(fs);
    return S;
}

void StructuresBinaryFile::jumpToStructureIndex(int idx) {
    MstUtils::assert(readMode, "jumpToStructureIndex not supported in write mode");
    if (idx < 0)
        fs.seekg(0, fs.beg);
    else if (idx >= _structureNames.size())
        fs.seekg(0, fs.end);
    else {
        string name = _structureNames[idx];
        fs.clear();
        fs.seekg(_filePositions[name], fs.beg);
    }
}

void StructuresBinaryFile::appendStructure(Structure *s) {
    MstUtils::assert(!readMode, "appendStructure not supported in read mode");
    s->writeData(fs);
}

// StructureCache

StructureCache::~StructureCache() {
    for (auto it = begin(); it != end(); ++it) {
        delete *it;
    }
}

void StructureCache::preloadFromBinaryFile() {
    MstUtils::assert(binaryFile != nullptr, "Cannot preload without a binary file");

    binaryFile->reset();
    while (binaryFile->hasNext() && (long)cache.size() < capacity) {
        Structure *s = binaryFile->next();
        cache.push_front(s);
        cachePointers[s->getName()] = cache.begin();
    }
    cout << "preload: load factor " << cachePointers.load_factor() << ", max " << cachePointers.max_load_factor() << endl;
}

Structure* StructureCache::getStructure(string name, string prefix) {
    Structure *structure = nullptr;
    auto pos = cachePointers.find(name);
    if (pos != cachePointers.end()) {
        structure = *(pos->second);

        // Remove the structure from the cache, since we'll put it at front
        cache.erase(pos->second);
    } else {
        if (binaryFile != nullptr) {
            structure = binaryFile->getStructureNamed(name);
        } else {
            string path = MstSystemExtension::join(prefix.size() > 0 ? prefix : pdbPrefix, name);
            structure = new Structure(path);
        }

        // Remove the back of the cache if it is full
        if ((long)cache.size() == capacity) {
            Structure *backStruct = cache.back();
            cache.erase(cachePointers[backStruct->getName()]);
            cachePointers.erase(backStruct->getName());
            atoms.erase(backStruct->getName());
            delete backStruct;
        }
    }

    // Push the new structure to the front
    cache.push_front(structure);
    cachePointers[name] = cache.begin();
    MstUtils::assert(structure->getName() == name, "Names don't match");
    return structure;
}

vector<Atom *> &StructureCache::getAtoms(string name, string prefix) {
    MstUtils::assert(storeAtoms, "storeAtoms must be true to use getAtoms");
    if (atoms.count(name) == 0) {
        Structure *s = getStructure(name, prefix);
        atoms[name] = s->getAtoms();
    }
    return atoms[name];
}

bool StructureCache::hasStructure(string name) {
    return cachePointers.count(name) != 0;
}

void StructureCache::removeStructure(string name) {
    if (cachePointers.count(name) != 0) {
        cache.erase(cachePointers[name]);
        cachePointers.erase(name);
    }
}

void StructureCache::clear() {
    cachePointers.clear();
    atoms.clear();
    cache.clear();
}

// StructureIterator

StructureIterator::StructureIterator(const vector<string>& filePaths, int batchSize, vector<string> *chainIDs): _filePaths(filePaths), _batchSize(batchSize), _chainIDs(chainIDs) {
    if (chainIDs != nullptr)
        MstUtils::assert(chainIDs->size() == filePaths.size(), "must have same number of chain IDs and file paths");
}

StructureIterator::StructureIterator(const string binaryFilePath, int batchSize, string chainID): _batchSize(batchSize), defaultChainID(chainID) {
    binaryFile = new StructuresBinaryFile(binaryFilePath);
}

bool StructureIterator::hasNext() {
    if (binaryFile != nullptr) {
        return binaryFile->hasNext();
    }
    return _batchIndex * _batchSize < _filePaths.size();
}

void StructureIterator::reset() {
    _batchIndex = 0;
    if (binaryFile != nullptr) {
        binaryFile->reset();
    }
}

vector<Structure *> StructureIterator::next() {
    int startIndex = (_batchIndex++) * _batchSize;
    vector<Structure *> result;
   
    if (binaryFile != nullptr) {
        while (binaryFile->hasNext() && result.size() < _batchSize) {
            Structure *candidate = binaryFile->next();
            Structure floatingChains;
            extractChains(*candidate, splitChainIDs(defaultChainID), floatingChains);
            result.push_back(new Structure(floatingChains)); 
            delete candidate;
        }
    } else {
        for (int i = startIndex; i < min(startIndex + _batchSize, (int)_filePaths.size()); i++) {
            string path = _filePaths[i];
            Structure candidate(path);
            Structure floatingChains;
            if (_chainIDs != nullptr) {
                extractChains(candidate, splitChainIDs((*_chainIDs)[i]), floatingChains);
            } else {
                floatingChains = candidate;
            }
            result.push_back(new Structure(floatingChains));
            /*if (result.size() % 2000 == 0)
                cout << result.size() << endl;*/
        }
    }
    
    // Remove last batch
    if (_lastBatch != nullptr)
        delete _lastBatch;
    _lastBatch = new vector<Structure *>(result);
    
    return result;
}

void StructureIterator::skip() {
    if (binaryFile != nullptr) {
        int i = 0;
        while (binaryFile->hasNext() && i < _batchSize) {
            binaryFile->skip();
            i++;
        }
        cout << "Skipped " << i << " structures" << endl;
    }
    _batchIndex++;
}

PairStructureIterator::PairStructureIterator(const vector<string>& filePaths, int batchSize, vector<string> *chainIDs): _first(filePaths, batchSize, chainIDs), _second(filePaths, batchSize, chainIDs) {}

PairStructureIterator::PairStructureIterator(const string binaryFilePath, int batchSize, string chainID): _first(binaryFilePath, batchSize, chainID), _second(binaryFilePath, batchSize, chainID) {}

bool PairStructureIterator::hasNext() {
    return _first.hasNext() || _second.hasNext();
}

void PairStructureIterator::reset() {
    _first.reset();
    _second.reset();
}

tuple<vector<Structure *>, vector<Structure *>> PairStructureIterator::next() {
    vector<Structure *> firstStructures;
    vector<Structure *> secondStructures;
    
    if (!_second.hasNext() || _currentFirst.size() == 0) {
        _second.reset();
        cout << "Loading from first" << endl;
        for (Structure *ptr: _currentFirst) {
            delete ptr;
        }
        firstStructures = _first.next();
    } else if (_currentFirst.size() > 0) {
        cout << "Using saved first structures" << endl;
        firstStructures = _currentFirst;
    }
    
    cout << "Loading from second, " << _second.hasNext() << endl;
    secondStructures = _second.next();
    _currentFirst = firstStructures;
    return make_tuple(firstStructures, secondStructures);
}

void PairStructureIterator::skip() {
    if (!_second.hasNext() || _currentFirst.size() == 0) {
        cout << "Skipping from first" << endl;
        _second.reset();
        for (Structure *ptr: _currentFirst) {
            delete ptr;
        }
        _currentFirst = _first.next();
    } else {
        cout << "Skipping from second" << endl;
        _second.skip();
    }
}

BatchPairStructureIterator::BatchPairStructureIterator(const string &binaryFilePath, int workerIndex, int numWorkers, int batchSize, string chainID): workerIndex(workerIndex), numWorkers(numWorkers), batchSize(batchSize), chainID(chainID) {
    binaryFile = new StructuresBinaryFile(binaryFilePath);
    binaryFile->scanFilePositions();
    binaryFile->reset();
    firstIndex = workerIndex - numWorkers;
    numRows = (int)ceil(binaryFile->structureCount() / batchSize);
    cout << "Work matrix has " << numRows << " rows" << endl;
}

BatchPairStructureIterator::~BatchPairStructureIterator() {
    if (binaryFile != nullptr)
        delete binaryFile;
}

bool BatchPairStructureIterator::hasNext() {
    if (!nextResultAvailable) {
        makeNextResult();
     }
     return nextResultAvailable;
}

void BatchPairStructureIterator::reset() {
    firstIndex = workerIndex - numWorkers;
    secondIndex = -1;
    for (Structure *s: currentFirst)
        delete s;
    for (Structure *s: currentSecond)
        delete s;
    currentFirst.clear();
    currentSecond.clear();
}

pair<vector<Structure *>, vector<Structure *>> BatchPairStructureIterator::next() {
    if (!nextResultAvailable) {
        makeNextResult();
    }
    nextResultAvailable = false;
    return make_pair(currentFirst, currentSecond);
}

void BatchPairStructureIterator::makeNextResult() {
    /*
     * We want to balance the work among numWorkers workers while 
     * minimizing file jumps and skips. The plan will be to divide the
     * square of tasks into rows, where each row corresponds to a batch of
     * first structures, and distribute the rows among the workers. Since
     * early rows will have fewer batches than later rows, distribute the
     * rows symmetrically as follows:
     *
     * 0|A
     * 1|BB
     * 2|CCC
     * 3|AAAA
     * 4|AAAAA
     * 5|CCCCCC
     * 6|BBBBBBB
     * 7|AAAAAAAA
     *  ---------
     *   01234567     
     */
    
    if (nextResultAvailable)
        return;

    // First make sure we're on the correct first index
    if (firstIndex <= secondIndex || !binaryFile->hasNext()) {
        if ((float)firstIndex < (float)numRows / 2.0f && (float)(firstIndex + numWorkers) >= (float)numRows / 2.0f)
            firstIndex = numRows - firstIndex - 1;
        else
            firstIndex += numWorkers;
        // Jump to the position in the file and load currentFirst
        cout << "Jumping to row " << firstIndex << endl;
        binaryFile->jumpToStructureIndex(firstIndex * batchSize);
        for (Structure *s: currentFirst)
            delete s;
        currentFirst.clear();
        while (binaryFile->hasNext() && currentFirst.size() < batchSize) {
            Structure *candidate = binaryFile->next();
            Structure floatingChains;
            extractChains(*candidate, splitChainIDs(chainID), floatingChains);
            currentFirst.push_back(new Structure(floatingChains)); 
            delete candidate;
        }

        if (currentFirst.empty())
            return;

        // Reset second index and move file handle to beginning
        secondIndex = -1;
        binaryFile->reset();
    }

    // Now add second structures
    secondIndex++;
    for (Structure *s: currentSecond)
        delete s;
    currentSecond.clear();
    while (binaryFile->hasNext() && currentSecond.size() < batchSize) {
        Structure *candidate = binaryFile->next();
        Structure floatingChains;
        extractChains(*candidate, splitChainIDs(chainID), floatingChains);
        currentSecond.push_back(new Structure(floatingChains)); 
        delete candidate;
    }

    nextResultAvailable = true;
}


