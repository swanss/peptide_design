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

bool StructuresBinaryFile::hasNext() {
    MstUtils::assert(readMode, "hasNext not supported in write mode");
    return fs.peek() != EOF;
}

Structure *StructuresBinaryFile::next() {
    MstUtils::assert(readMode, "next not supported in write mode");
    return readNextFileSection(false).first;
}

void StructuresBinaryFile::skip() {
    MstUtils::assert(readMode, "skip not supported in write mode");
    Structure* S = readNextFileSection(false).first;
    delete S;
}

void StructuresBinaryFile::reset() {
    MstUtils::assert(readMode, "reset not supported in write mode");
    fs.clear(); // this is necessary in case ifstream doesn't clear eofbit
    fs.seekg(0, fs.beg);
}

void StructuresBinaryFile::scanFilePositions() {
    //MstUtils::assert(readMode, "scanFilePositions not supported in write mode");
    MstTimer timer; timer.start();
    if (!_structureNames.empty()) return;
    reset();
    cout << "Scanning file positions..." << endl;
    _structureNames.clear();
    int count = 0;
    while (fs.peek() != EOF) {
        pair<Structure*,long> next;
        next = readNextFileSection(true,false);
    //    if (count % 10000 == 0) next = readNextFileSection(true,true);
    //    else next = readNextFileSection(true);
        Structure *S = next.first;
        long pos = next.second;
        if (pos < 0) {
            cout << pos << endl;
        }
        _filePositions[S->getName()] = pos;
        _structureNames.push_back(S->getName());
        // if (count % 10000 == 0) {
        //     // report values
        //     cout << "Structure number: " << count << endl;
        //     cout << "Number of structure names:" << _structureNames.size() << endl;
        //     cout << "Number of file positions: " << _filePositions.size() << endl;
        //     cout << "Number of seed int val types: " << seed_dscrt_vals.size() << endl;
        //     for (auto it : seed_dscrt_vals) {
        //         cout << it.first << " with size: " << it.second.size() << endl;
        //     }
        //     cout << "Number of seed real val types: " << seed_real_vals.size() << endl;
        //     for (auto it : seed_real_vals) {
        //         cout << it.first << " with size: " << it.second.size() << endl;
        //     }
        // }
        delete S;
        count++;
    }
    timer.stop();
    cout << "Done scanning file, took " << timer.getDuration(MstTimer::msec) / 1000.0 << " sec" << endl;
}

Structure * StructuresBinaryFile::getStructureNamed(string name) {
    MstUtils::assert(readMode, "getStructureNamed not supported in write mode");
    if (_filePositions.size() == 0) {
        scanFilePositions();
    }
    if (_filePositions.count(name) == 0) {
        cout << "Structure " << name << " doesn't exist" << endl;
        return nullptr;
    }
    fs.seekg(_filePositions[name], fs.beg);
    return readNextFileSection(false).first;
}

void StructuresBinaryFile::jumpToStructureIndex(int idx) {
    MstUtils::assert(readMode, "jumpToStructureIndex not supported in write mode");
    if (_structureNames.size() == 0) {
        scanFilePositions();
    }
    if (idx < 0) {
        fs.seekg(0, fs.beg);
    } else if (idx >= _structureNames.size()) {
        fs.seekg(0, fs.end);
    } else {
        string name = _structureNames[idx];
        fs.clear();
        fs.seekg(_filePositions[name], fs.beg);
    }
}

void StructuresBinaryFile::appendStructure(Structure *s) {
    MstUtils::assert(!readMode, "appendStructure not supported in read mode");
    MstUtils::assert(s->residueSize() != 0, "Structure must have at least one residue");

    if (!structure_added) {
        structure_added = true;
        if (!_append) MstUtils::writeBin(fs,string("!@#version_2!@#")); //write the version at the top of the file
    } else {
        MstUtils::writeBin(fs,'E'); //finish the previous section
    }
    MstUtils::writeBin(fs,'S'); //start new structure section
    s->writeData(fs);
}

void StructuresBinaryFile::appendStructurePropertyInt(string prop, int val) {
    MstUtils::assert(!prop.empty(), "property field cannot be empty");

    MstUtils::writeBin(fs,'I');
    MstUtils::writeBin(fs,prop);
    MstUtils::writeBin(fs,val);
}

void StructuresBinaryFile::appendStructurePropertyReal(string prop, mstreal val) {
    MstUtils::assert(!prop.empty(), "property field cannot be empty");
    
    MstUtils::writeBin(fs,'R');
    MstUtils::writeBin(fs,prop);
    MstUtils::writeBin(fs,val);
}

void StructuresBinaryFile::openFileStream() {
    if (readMode)
        MstUtils::openFile(fs, _filePath, fstream::in | fstream::binary, "StructuresBinaryFile::openFileStream");
    else if (_append)
        MstUtils::openFile(fs, _filePath, fstream::out | fstream::binary | fstream::app, "StructuresBinaryFile::openFileStream");
    else MstUtils::openFile(fs, _filePath, fstream::out | fstream::binary, "StructuresBinaryFile::openFileStream");
}

void StructuresBinaryFile::detectFileVersion() {
    if (readMode == false) {
        _version = 2;
        return;
    }
    MstUtils::openFile(fs, _filePath, fstream::in | fstream::binary, "StructuresBinaryFile::openFileStream");
    string version;
    MstUtils::readBin(fs, version);
    if (version == "!@#version_2!@#") _version = 2;
    else _version = 1;
    fs.close();
}

pair<Structure*,long> StructuresBinaryFile::readNextFileSection(bool save_metadata, bool verbose) {
    //if beginning of file, advance past the version
    if ((_version == 2) && (fs.tellg() == 0)) {
        string version; MstUtils::readBin(fs, version);
    }
    Structure* S = new Structure();
    long pos = fs.tellg();
    if (_version == 1) {
        S->readData(fs);
        return pair<Structure*,long>(S,pos);
    } else if (_version == 2) {
        char sect; string prop; mstreal real_val; int dscrt_val;
        MstUtils::readBin(fs, sect);
        if (sect != 'S') MstUtils::error("The first section should be a Structure " + _filePath, "StructuresBinaryFile::readFileSection()");
        S->readData(fs);
        if (verbose) cout << S->getName() << endl;
        //read meta-data;
        while (fs.peek() != EOF) {
            MstUtils::readBin(fs, sect);
            if (sect == 'I') {
                MstUtils::readBin(fs, prop);
                MstUtils::readBin(fs, dscrt_val);
                if (save_metadata) seed_dscrt_vals[prop][S->getName()] = dscrt_val;
                if (verbose) cout << prop << " : " << dscrt_val << endl;
            } else if (sect == 'R') {
                MstUtils::readBin(fs,prop);
                MstUtils::readBin(fs, real_val);
                if (save_metadata) seed_real_vals[prop][S->getName()] = real_val;
                if (verbose) cout << prop << " : " << real_val << endl;
            } else if (sect == 'E') {
                break;
            } else {
                MstUtils::error("Section name not recognized: " + MstUtils::toString(sect),"StructuresBinaryFile::readFileSection()");
            }
        }
        return pair<Structure*,long>(S,pos);
    }
    return pair<Structure*,long>(S,pos); //just to quell compiler warning;
}

vector<string> StructuresBinaryFile::getStructureNames() {
    MstUtils::assert(readMode, "getStructureNamed not supported in write mode");
    if (_filePositions.size() == 0) {
        scanFilePositions();
    }
    return _structureNames;
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
    preloaded = true;
    if (cache.size() < capacity) _belowCapacity = true;
}

void StructureCache::preloadFromPDBList(string pdbList) {
    vector<string> pdbNames = MstUtils::fileToArray(pdbList);
    for (string pdbName : pdbNames) {
        string path = MstSystemExtension::join(pdbPrefix, pdbName);
        Structure* s = new Structure(path);
        cache.push_front(s);
        cachePointers[s->getName()] = cache.begin();
        if ((long)cache.size() == capacity) break;
    }
    cout << "preload: load factor " << cachePointers.load_factor() << ", max " << cachePointers.max_load_factor() << endl;
    preloaded = true;
    if (cache.size() < capacity) _belowCapacity = true;
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
    MstUtils::assert(structure->getName() == name, "Names don't match: "+structure->getName()+" vs "+name);
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

StructureIterator::StructureIterator(const vector<string>& filePaths, int batchSize, vector<string> *chainIDs, int workerIndex, int numWorkers): _filePaths(filePaths), _batchSize(batchSize), _chainIDs(chainIDs), _workerIndex(workerIndex), _numWorkers(numWorkers), _batchIndex(workerIndex) {
    if (chainIDs != nullptr)
        MstUtils::assert(chainIDs->size() == filePaths.size(), "must have same number of chain IDs and file paths");
}

StructureIterator::StructureIterator(const string binaryFilePath, int batchSize, string chainID, int workerIndex, int numWorkers): _batchSize(batchSize), defaultChainID(chainID), _workerIndex(workerIndex), _numWorkers(numWorkers), _batchIndex(workerIndex) {
    binaryFile = new StructuresBinaryFile(binaryFilePath);
}

bool StructureIterator::hasNext() {
    if (!nextResultAvailable) makeNextResult();
    return nextResultAvailable;
}

void StructureIterator::reset() {
    _batchIndex = _workerIndex;
    if (binaryFile != nullptr) {
        binaryFile->reset();
    }
}

vector<Structure *> StructureIterator::next() {
    if (!nextResultAvailable) makeNextResult();
    nextResultAvailable = false;
    return *_currentBatch;
}

/**
 1) Find the structure index using the batch index
 2) Load the structures
 2) Advance the batch index
 */

void StructureIterator::makeNextResult() {
    vector<Structure *> result;
    
    // get structure index using the batch index
    int startIndex = (_batchIndex - 1) * _batchSize;
    
    if (binaryFile != nullptr) {
        // go to start index in binary file
        binaryFile->jumpToStructureIndex(startIndex);
        
        while (binaryFile->hasNext() && result.size() < _batchSize) {
            Structure *candidate = binaryFile->next();
            string name = candidate->getName();
            Structure floatingChains;
            extractChains(*candidate, splitChainIDs(defaultChainID), floatingChains);
            Structure* newStructure = new Structure(floatingChains);
            newStructure->setName(name);
            result.push_back(newStructure);
            
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
            
            Structure* newStructure = new Structure(floatingChains);
            newStructure->setName(path);
            result.push_back(newStructure);
            /*if (result.size() % 2000 == 0)
                cout << result.size() << endl;*/
        }
    }
    
    // Replace the current batch with the new structures
    if (_currentBatch != nullptr) {
        if (maintainsOwnership) for (Structure* seed : *_currentBatch) delete seed;
        delete _currentBatch;
    }
    _currentBatch = new vector<Structure *>(result);
    if (!_currentBatch->empty()) nextResultAvailable = true;
    
    // Advance the batch index
    _batchIndex = _batchIndex + _numWorkers;
}

void StructureIterator::skip() {
    if (binaryFile != nullptr) {
        int startIndex = (_batchIndex - 1) * _batchSize;
        binaryFile->jumpToStructureIndex(startIndex);
    }
    
    //increment the batch index
    _batchIndex = _batchIndex + _numWorkers;
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
    numRows = (int)ceil((double)binaryFile->structureCount() / (double)batchSize);
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
    // Ensure that if the first and second parts of the batch are the same structures,
    // the same underlying objects are returned
    if (currentSecond[0]->getName() == currentFirst[0]->getName() && currentSecond.size() == currentFirst.size())
        return make_pair(currentFirst, currentFirst);
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
        // this condition checks if firstIndex would cross the halfway point (numRows / 2)
        // in this iteration of the loop (only applicable if there aren't enough workers to
        // cover all rows)
        if (numWorkers < numRows / 2 && firstIndex < numRows / 2 && (firstIndex + numWorkers) >= numRows / 2) {
            // If the halfway point will be crossed, the order of assigning batches is flipped.
            // The next line looks at the distance to the halfway point and treats it like a mirror.
            firstIndex = firstIndex + 2 * (numRows / 2 - firstIndex) - 1;
        }
        else firstIndex += numWorkers;
        
        // Jump to the position in the file and load currentFirst
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

        cout << "Starting work row " << firstIndex + 1 << "/" << numRows << endl;
        
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


