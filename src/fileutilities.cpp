//
//  fileformats.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 1/28/19.
//

#include <stdio.h>
#include "fileutilities.h"
#include "utilities.h"

// FuseCandidateFile

FuseCandidateFile::~FuseCandidateFile() {
    if (readstream != nullptr) {
        delete readstream;
    }
    if (writestream != nullptr) {
        delete writestream;
    }
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
