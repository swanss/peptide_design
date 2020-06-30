// Implementation for OverlapFinder template class
#include <algorithm>

template <class Hasher>
void OverlapFinder<Hasher>::insertStructures(const vector<Structure *> &structures) {
    // Hash every residue that is in a chain with the correct ID
    for (Structure *s: structures) {
        for (int chainIdx = 0; chainIdx < s->chainSize(); ++chainIdx) {
            Chain &c = s->getChain(chainIdx);
            if (!_seedChain.empty() && c.getID() != _seedChain)
                continue;
            
            auto residues = c.getResidues();
            for (int resIdx = 0; resIdx < residues.size(); ++resIdx) {
                auto hashValue = _hasher.hash(residues[resIdx]);
                _hashMap[hashValue].emplace_back(residues[resIdx], resIdx);
            }
        }
    }
}

template <class Hasher>
vector<FuseCandidate> OverlapFinder<Hasher>::findOverlaps(const vector<Structure *> &testStructures, bool constrainOrder) const {
    
    vector<FuseCandidate> results;
    
    // Loop over appropriate chains in all structures to test
    for (int i = 0; i < testStructures.size(); i++) {
        if (i % (testStructures.size() / 100) == 0)
            cout << (i / (testStructures.size() / 100)) << "% complete (" << i << "/" << testStructures.size() << ") - " << results.size() << " overlaps" << endl;
        Structure *s = testStructures[i];
        for (int chainIdx = 0; chainIdx < s->chainSize(); ++chainIdx) {
            Chain &c = s->getChain(chainIdx);
            if (!_seedChain.empty() && c.getID() != _seedChain)
                continue;
            
            auto residues = c.getResidues();
            findOverlaps(residues, results, constrainOrder);
        }
    }
    
    return results;
}

template <class Hasher>
void OverlapFinder<Hasher>::findOverlaps(const vector<Structure *> &testStructures, FuseCandidateFile &outFile, bool constrainOrder) const {
    
    int numResults = 0;
    
    // Loop over appropriate chains in all structures to test
    for (int i = 0; i < testStructures.size(); i++) {
        if (i % (testStructures.size() / 100) == 0)
            cout << (i / (testStructures.size() / 100)) << "% complete (" << i << "/" << testStructures.size() << ") - " << numResults << " overlaps" << endl;
        Structure *s = testStructures[i];
        
        vector<FuseCandidate> results;
        for (int chainIdx = 0; chainIdx < s->chainSize(); ++chainIdx) {
            Chain &c = s->getChain(chainIdx);
            if (!_seedChain.empty() && c.getID() != _seedChain)
                continue;
            
            auto residues = c.getResidues();
            findOverlaps(residues, results, constrainOrder);
        }
        numResults += results.size();
        
        // Write results to file
        outFile.write(results, "");
    }
}

template <class Hasher>
void OverlapFinder<Hasher>::findOverlaps(vector<Residue *> &source, vector<FuseCandidate> &results, bool constrainOrder) const {
    
    // Add overlap segment candidates for each residue in the segment
    vector<vector<OverlapCandidate>> candidates;
    // Keep a list of pointers to the beginning of each vector (for later)
    typedef vector<OverlapCandidate>::const_iterator Cursor;
    vector<Cursor> cursors;
    
    int segmentIndex = 0;
    for (Residue *res: source) {
        candidates.emplace_back();
        vector<OverlapCandidate> &positionCandidates = candidates.back();
        
        Structure *parentStructure = res->getStructure();
        
        vector<typename Hasher::hash_type> region = _hasher.region(res, _cutoff);
        for (typename Hasher::hash_type &hashValue: region) {
            // Add overlap segments for all the residues in each region
            if (_hashMap.count(hashValue) == 0)
                continue;
            
            for (const pair<Residue *, int> &candidate: _hashMap.at(hashValue)) {
                // Exclude residues from the same structure. If we're constraining the
                // order so duplicate overlaps aren't returned, simply check the order
                // of the structures' pointer addresses. Since we assume all
                // structures will use the same addresses throughout the program, this
                // is an unambiguous ordering.
                Residue *candRes = candidate.first;
                Structure *candStructure = candRes->getStructure();
                if (candStructure == parentStructure || (constrainOrder && candStructure < parentStructure))
                    continue;

                positionCandidates.emplace_back(candStructure, candRes->getParent()->getID(), candidate.second);
            }
            if (positionCandidates.size() >= 1000)
                break;
        }
        
        // The algorithm depends on the candidates for each position being sorted
        // in a consistent way.
        sort(positionCandidates.begin(), positionCandidates.end(), OverlapCandidate::compare);
        
        cursors.push_back(positionCandidates.begin());
        
        segmentIndex++;
    }
    
    // Can't do anything if there are no residues in the segment
    if (candidates.empty())
        return;
    
    // Now, iterate through all lists of candidates and look for occurrences of
    // the same overlap segment in a contiguous sequence.
    while (true) {
        // The only overlap we consider in each iteration is the segment with the "min"
        // value, since it's guaranteed that we won't see that segment ever again.
        const OverlapCandidate *minElement = nullptr;
        for (int i = 0; i < candidates.size(); i++) {
            if (cursors[i] == candidates[i].end()) continue;
            if (minElement == nullptr || OverlapCandidate::compare(*cursors[i], *minElement)) {
                minElement = &(*cursors[i]);
            }
            
        }
        
        vector<Cursor> overlapPositions;
        // Keep track of whether there are still contiguous segments with available overlaps
        int contiguousSection = 0;
        bool hasContiguousSegment = false;
        
        // Look for sets of positions that match minElement with consecutive residue indexes
        for (int i = 0; i < candidates.size(); i++) {
            if (cursors[i] == candidates[i].end()) {
                contiguousSection = 0;
                continue;
            }
            contiguousSection++;
            if (contiguousSection >= _segmentLength)
                hasContiguousSegment = true;
            
            // If the current position doesn't match the min element, reset and skip this position
            if (!cursors[i]->fromSameChain(*minElement)) {
                overlapPositions.clear();
                continue;
            }
            
            // If the current position is noncontiguous, ignore it
            if (!overlapPositions.empty() && cursors[i]->offset != overlapPositions.back()->offset + 1) {
                overlapPositions.clear();
            }
            
            overlapPositions.push_back(cursors[i]);
            if (overlapPositions.size() >= _segmentLength) {
                // A possible overlap!
                const OverlapCandidate &startPosition = *overlapPositions[overlapPositions.size() - _segmentLength];

                if (_verifier == nullptr || checkOverlap(source.begin() + i + 1 - _segmentLength, source.begin() + i + 1, startPosition)) {
                    // Build a fuse candidate object
                    FuseCandidate fuseCandidate;
                    fuseCandidate.overlapSize = _segmentLength;
                    Residue *firstRes = source[i + 1 - _segmentLength];
                    fuseCandidate.setStructure1(firstRes->getStructure(), firstRes->getParent()->getID());
                    fuseCandidate.overlapPosition1 = i + 1 - _segmentLength;
                    fuseCandidate.setStructure2(startPosition.structure, startPosition.chainID);
                    fuseCandidate.overlapPosition2 = startPosition.offset;
                    
                    // Don't store RMSD in this implementation
                    fuseCandidate.rmsd = 0.0;
                    results.push_back(fuseCandidate);
                }
            }
        }
        
        // If there is no stretch of _segmentLength positions that still has candidates available,
        // there can be no more overlaps
        if (!hasContiguousSegment)
            break;
        
        bool foundCursor = false;
        for (int i = 0; i < candidates.size(); i++) {
            if (cursors[i] == candidates[i].end()) continue;
            if ((*cursors[i]).fromSameChain(*minElement)) {
                cursors[i]++;
            }
            if (cursors[i] != candidates[i].end())
                foundCursor = true;
        }
        
        // If all cursors have reached the end, we're done
        if (!foundCursor)
            break;
    }
}

template <class Hasher>
bool OverlapFinder<Hasher>::checkOverlap(vector<Residue *>::iterator begin, vector<Residue *>::iterator end, const OverlapCandidate &candidate) const {
    if (!_verifier)
        return true;
                
    // Build lists of residues for each segment
    vector<Residue *> segment1(begin, end);
    
    vector<Residue *> chainResidues = candidate.structure->getChainByID(candidate.chainID)->getResidues();
    vector<Residue *> segment2(chainResidues.begin() + candidate.offset, chainResidues.begin() + candidate.offset + _segmentLength);
    
    // Verify
    return _verifier->verify(segment1, segment2);
}
