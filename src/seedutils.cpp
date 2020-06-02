//
//  fileformats.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 1/28/19.
//

#include <stdio.h>
#include "seedutils.h"
#include <regex>

int getTargetResidueIndex(string seedName) {
    regex re("^[A-z0-9]{4}_[A-z](\\d+)_");
    smatch match;
    string fileName = MstSystemExtension::fileName(seedName);
    if (std::regex_search(fileName, match, re) && match.size() > 1) {
        return atoi(match.str(1).c_str());
    }
    return -1;
}

pair<string, int> getTargetResidueCode(string seedName) {
    regex re("^[A-z0-9_]-([A-z])(\\d+)-");
    smatch match;
    string fileName = MstSystemExtension::fileName(seedName);
    if (std::regex_search(fileName, match, re) && match.size() > 1) {
        return make_pair(match.str(1), atoi(match.str(2).c_str()));
    }
    return make_pair("", -1);
}

// FuseCandidateFile

FuseCandidateFile::~FuseCandidateFile() {
    if (readstream != nullptr) {
        delete readstream;
    }
    if (writestream != nullptr) {
        delete writestream;
    }
}

vector<string> splitString(string s, const string &delim) {
    size_t pos = 0;
    string token;
    vector<string> result;
    while ((pos = s.find(delim)) != string::npos) {
        token = s.substr(0, pos);
        result.push_back(token);
        s.erase(0, pos + delim.length());
    }
    result.push_back(s);
    return result;
}

vector<FuseCandidate> FuseCandidateFile::read(int numLines) {
    vector<FuseCandidate> results;
    
    if (readstream == nullptr) {
        if (writestream != nullptr) {
            cerr << "Can't read and write from same file" << endl;
        }
        readstream = new ifstream(_path);
        if (!readstream->is_open()) {
            cerr << "couldn't open in stream" << endl;
            return results;
        }
    }
    
    string line;
    while (getline(*readstream, line)) {
        vector<string> comps = splitString(line, ",");
        FuseCandidate candidate;
        candidate.file1 = comps[0];
        candidate.chain1 = comps[1];
        candidate.file2 = comps[2];
        candidate.chain2 = comps[3];
        candidate.overlapPosition1 = stoi(comps[4]);
        candidate.overlapPosition2 = stoi(comps[5]);
        candidate.overlapSize = stoi(comps[6]);
        candidate.rmsd = stof(comps[7]);
        results.push_back(candidate);
        if (results.size() == numLines)
            break;
    }
    return results;
}

string getRelativePathIfPossible(string path, string parent) {
    if (parent.size() == 0)
        return MstSystemExtension::fileName(path);
    return MstSystemExtension::relativePath(path, parent);
}

void FuseCandidateFile::write(vector<FuseCandidate> candidates, string parentPath) {
    if (writestream == nullptr) {
        if (readstream != nullptr) {
            cerr << "Can't read and write from same file" << endl;
        }
        writestream = new ofstream(_path, _shouldAppend ? (ofstream::app | ofstream::out) : ios::out);
        if (!writestream->is_open()) {
            cerr << "couldn't open out stream" << endl;
            return;
        }
    }
    
    for (FuseCandidate candidate: candidates) {
        *writestream << getRelativePathIfPossible(candidate.file1, parentPath) << ","
        << candidate.chain1 << ","
        << getRelativePathIfPossible(candidate.file2, parentPath) << ","
        << candidate.chain2 << ","
        << candidate.overlapPosition1 << ","
        << candidate.overlapPosition2 << ","
        << candidate.overlapSize << ","
        << candidate.rmsd << endl;
    }
}

void FuseCandidateFile::write(FuseCandidate candidate, string parentPath) {
    vector<FuseCandidate> candidates;
    candidates.push_back(candidate);
    write(candidates, parentPath);
}

// SeedListFile

SeedListFile::~SeedListFile() {
    if (readstream != nullptr) {
        delete readstream;
    }
    if (writestream != nullptr) {
        delete writestream;
    }
}

pair<vector<string>, vector<string>> SeedListFile::read(string pdbPrefix) {
    vector<string> results;
    vector<string> chainIDs;
    
    if (readstream == nullptr) {
        if (writestream != nullptr) {
            cerr << "Can't read and write from same file" << endl;
        }
        readstream = new ifstream(_path);
        if (!readstream->is_open()) {
            cerr << "couldn't open in stream" << endl;
            return pair<vector<string>, vector<string>>(results, chainIDs);
        }
    }
    
    string line;
    while (getline(*readstream, line)) {
        vector<string> comps = splitString(line, ",");
        string path = comps[0];
        if (path == "path") {
            // Header row
            MstUtils::assert(comps[1] == "chids", "expected second column to be 'chids'");
            _metadataNames.clear();
            _metadataNames.insert(_metadataNames.end(), comps.begin() + 2, comps.end());
        } else {
            for (int i = 2; i < comps.size(); i++) {
                metadata[path].push_back(comps[i]);
            }
            results.push_back(pdbPrefix + path);
            chainIDs.push_back(comps[1]);
        }
    }
    return pair<vector<string>, vector<string>>(results, chainIDs);
}

void SeedListFile::write(string seed, string chainIDs, vector<string> metadata) {
    if (writestream == nullptr) {
        if (readstream != nullptr) {
            cerr << "Can't read and write from same file" << endl;
        }
        writestream = new ofstream(_path, ios::out);
        if (!writestream->is_open()) {
            cerr << "couldn't open out stream" << endl;
            return;
        }
        
        // Write header
        *writestream << "path,chids";
        for (string headerName: _metadataNames) {
            *writestream << "," << headerName;
        }
        *writestream << endl;
    }
    
    *writestream << seed << "," << chainIDs;
    MstUtils::assert(metadata.size() == _metadataNames.size(), metadata.size() == 0 ? "must provide metadata if field names are provided!" : "must provide metadata field names before writing a seed with metadata!");
    for (string field: metadata) {
        *writestream << "," << field;
    }
    *writestream << endl;
}

template <class T>
long PositionMap<T>::hash(float x, float y, float z) {
    MstUtils::assert(fabs(x) < MAX_COORDINATE && fabs(y) < MAX_COORDINATE && fabs(z) < MAX_COORDINATE, "coordinates out of bounds");
    
    long xBin = (long)floor(x / _interval) + _maxCoordBin / 2;
    long yBin = (long)floor(y / _interval) + _maxCoordBin / 2;
    long zBin = (long)floor(z / _interval) + _maxCoordBin / 2;
    cout << xBin << ", " << yBin << ", " << zBin << endl;
    return xBin * _maxCoordBin * _maxCoordBin + yBin * _maxCoordBin + zBin;
}

template <class T>
vector<long> PositionMap<T>::hash(float xmin, float ymin, float zmin, float xmax, float ymax, float zmax) {
    MstUtils::assert(fabs(xmin) < MAX_COORDINATE && fabs(ymin) < MAX_COORDINATE && fabs(zmin) < MAX_COORDINATE && fabs(xmax) < MAX_COORDINATE && fabs(ymax) < MAX_COORDINATE && fabs(zmax) < MAX_COORDINATE, "coordinates out of bounds");
    
    long xBinMin = (long)floor(xmin / _interval) + _maxCoordBin / 2;
    long yBinMin = (long)floor(ymin / _interval) + _maxCoordBin / 2;
    long zBinMin = (long)floor(zmin / _interval) + _maxCoordBin / 2;
    long xBinMax = (long)ceil(xmax / _interval) + _maxCoordBin / 2;
    long yBinMax = (long)ceil(ymax / _interval) + _maxCoordBin / 2;
    long zBinMax = (long)ceil(zmax / _interval) + _maxCoordBin / 2;
    vector<long> result;
    for (long x = xBinMin; x < xBinMax; x++) {
        for (long y = yBinMin; y < yBinMax; y++) {
            for (long z = zBinMin; z < zBinMax; z++) {
                result.push_back(x * _maxCoordBin * _maxCoordBin + y * _maxCoordBin + z);
            }
        }
    }
    return result;
}

template <class T>
void PositionMap<T>::insert(float x, float y, float z, T identifier) {
    _table[hash(x, y, z)].push_back(identifier);
}

template <class T>
void PositionMap<T>::insert(float xmin, float ymin, float zmin, float xmax, float ymax, float zmax, T identifier) {
    for (long bin: hash(xmin, ymin, zmin, xmax, ymax, zmax)) {
        _table[bin].push_back(identifier);
    }
}

template <class T>
vector<T> PositionMap<T>::itemsNear(float x, float y, float z) {
    vector<T> result;
    for (long bin: nearbyBins(x, y, z)) {
        auto it = _table.find(bin);
        if (it != _table.end()) {
            vector<T> nearby = it->second;
            result.insert(result.end(), nearby.begin(), nearby.end());
        }
    }
    return result;
}

template <class T>
void PositionMap<T>::pushIfInRange(long bin, vector<long> &arr) {
    if (bin >= 0 && bin < _maxBin)
        arr.push_back(bin);
}

template <class T>
vector<long> PositionMap<T>::nearbyBins(float x, float y, float z) {
    long bin = hash(x, y, z);
    vector<long> result;
    pushIfInRange(bin, result);
    pushIfInRange(bin + 1, result);
    pushIfInRange(bin - 1, result);
    pushIfInRange(bin + _maxCoordBin, result);
    pushIfInRange(bin - _maxCoordBin, result);
    pushIfInRange(bin + _maxCoordBin + 1, result);
    pushIfInRange(bin + _maxCoordBin - 1, result);
    pushIfInRange(bin - _maxCoordBin + 1, result);
    pushIfInRange(bin - _maxCoordBin - 1, result);
    
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin, result);
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin + 1, result);
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin - 1, result);
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin + _maxCoordBin, result);
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin - _maxCoordBin, result);
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin + _maxCoordBin + 1, result);
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin + _maxCoordBin - 1, result);
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin - _maxCoordBin + 1, result);
    pushIfInRange(bin + _maxCoordBin * _maxCoordBin - _maxCoordBin - 1, result);
    
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin, result);
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin + 1, result);
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin - 1, result);
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin + _maxCoordBin, result);
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin - _maxCoordBin, result);
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin + _maxCoordBin + 1, result);
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin + _maxCoordBin - 1, result);
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin - _maxCoordBin + 1, result);
    pushIfInRange(bin - _maxCoordBin * _maxCoordBin - _maxCoordBin - 1, result);
    return result;
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
            cout << "Loading structure" << endl;
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
