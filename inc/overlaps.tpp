// Implementation for OverlapFinder template class
#include <algorithm>

template <class Hasher>
void OverlapFinder<Hasher>::insertStructures(const vector<Structure *> &structures) {
    // Hash every residue that is in a chain with the correct ID
    cout << "Inserting " << structures.size() << " structures" << endl;
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
    cout << "Done inserting" << endl;
}

template <class Hasher>
vector<FuseCandidate> OverlapFinder<Hasher>::findOverlaps(const vector<Structure *> &testStructures) const {
    
    vector<FuseCandidate> results;
    
    // Loop over appropriate chains in all structures to test
    for (int i = 0; i < testStructures.size(); i++) {
        if (i % 100 == 0)
            cout << "Structure " << i << " of " << testStructures.size() << endl;
        Structure *s = testStructures[i];
        for (int chainIdx = 0; chainIdx < s->chainSize(); ++chainIdx) {
            Chain &c = s->getChain(chainIdx);
            if (!_seedChain.empty() && c.getID() != _seedChain)
                continue;
            
            auto residues = c.getResidues();
            
            // Iterate over all segments of length _segmentLength
            int pos = 0;
            for (auto it = residues.begin(); it != residues.begin() + residues.size() - _segmentLength + 1; ++it) {
                findOverlaps(it, it + _segmentLength, pos++, results);
            }
        }
    }
    
    return results;
}

template <class Hasher>
void OverlapFinder<Hasher>::findOverlaps(const vector<Structure *> &testStructures, FuseCandidateFile &outFile) const {
    
    // Loop over appropriate chains in all structures to test
    for (int i = 0; i < testStructures.size(); i++) {
        if (i % 100 == 0)
            cout << "Structure " << i << " of " << testStructures.size() << endl;
        Structure *s = testStructures[i];
        for (int chainIdx = 0; chainIdx < s->chainSize(); ++chainIdx) {
            Chain &c = s->getChain(chainIdx);
            if (!_seedChain.empty() && c.getID() != _seedChain)
                continue;
            
            auto residues = c.getResidues();
            
            // Iterate over all segments of length _segmentLength
            int pos = 0;
            vector<FuseCandidate> results;
            for (auto it = residues.begin(); it != residues.begin() + residues.size() - _segmentLength + 1; ++it) {
                findOverlaps(it, it + _segmentLength, pos++, results);
            }
            
            // Write results to file
            outFile.write(results, "");
        }
    }
}

template <class Hasher>
void OverlapFinder<Hasher>::findOverlaps(vector<Residue *>::iterator begin, vector<Residue *>::iterator end, int firstOffset, vector<FuseCandidate> &results) const {
    
    // Add overlap segment candidates for each residue in the segment
    vector<vector<OverlapSegment>> candidates;
    // Keep a list of pointers to the beginning of each vector (for later)
    typedef vector<OverlapSegment>::const_iterator Cursor;
    vector<Cursor> cursors;
    
    int segmentIndex = 0;
    for (auto it = begin; it != end; ++it) {
        candidates.emplace_back();
        vector<OverlapSegment> &positionCandidates = candidates.back();
        
        Residue *res = *it;
        Structure *parentStructure = res->getStructure();
        
        vector<typename Hasher::hash_type> region = _hasher.region(res, _cutoff);
        for (typename Hasher::hash_type &hashValue: region) {
            // Add overlap segments for all the residues in each region
            if (_hashMap.count(hashValue) == 0)
                continue;
            
            for (const pair<Residue *, int> &candidate: _hashMap.at(hashValue)) {
                // Exclude residues from the same structure, and from structures that
                // have pointer addresses less than the reference. Since we assume all
                // structures will use the same addresses throughout the program, this
                // is an unambiguous ordering.
                Residue *candRes = candidate.first;
                Structure *candStructure = candRes->getStructure();
                if (candStructure <= parentStructure)
                    continue;
                
                int offset = candidate.second - segmentIndex;
                // If the offset needed for this overlap is negative, we won't be
                // able to form a long enough overlap anyway
                if (offset < 0)
                    continue;
                
                positionCandidates.emplace_back(candStructure, candRes->getParent()->getID(), offset);
            }
        }
        
        // If there are no candidates for any one position, there can be no overlaps
        if (positionCandidates.empty())
            return;
        
        // The algorithm depends on the candidates for each position being sorted
        // in a consistent way.
        sort(positionCandidates.begin(), positionCandidates.end(), OverlapSegment::compare);
        
        cursors.push_back(positionCandidates.begin());
        
        segmentIndex++;
    }
    
    // Can't do anything if there are no residues in the segment
    if (candidates.empty())
        return;
    
    // Now, iterate jointly through all lists of candidates and look for occurrences of
    // the same overlap segment in *all* candidate lists. This is akin to a set
    // intersection operation.
    while (true) {
        // Compare all current cursor positions to the first one (they all need to be
        // equal)
        const OverlapSegment &referenceItem = *cursors[0];
        bool matched = true;
        for (const Cursor &cursor: cursors) {
            if (*cursor != referenceItem) {
                matched = false;
                break;
            }
        }
        
        if (matched) {
            // Great - a potential overlap! Verify it and add it to the list
            
            if (_verifier == nullptr || checkOverlap(begin, end, referenceItem)) {
                // Build a fuse candidate object
                FuseCandidate fuseCandidate;
                fuseCandidate.overlapSize = _segmentLength;
                Residue *firstRes = *begin;
                fuseCandidate.setStructure1(firstRes->getStructure(), firstRes->getParent()->getID());
                fuseCandidate.overlapPosition1 = firstOffset;
                fuseCandidate.setStructure2(referenceItem.structure, referenceItem.chainID);
                fuseCandidate.overlapPosition2 = referenceItem.offset;
                
                // Don't store RMSD in this implementation
                fuseCandidate.rmsd = 0.0;
                results.push_back(fuseCandidate);
            }
            
            // Advance all of the iterators
            bool terminate = false;
            for (int i = 0; i < cursors.size(); i++) {
                cursors[i]++;
                
                // Stop the whole loop when any one iterator reaches the end
                if (cursors[i] == candidates[i].end()) {
                    terminate = true;
                    break;
                }
            }
            
            if (terminate)
                break;
            
        } else {
            // Not an overlap - advance the cursor that points to the minimum sorted value
            int minIndex = min_element(cursors.begin(), cursors.end(), +[](const Cursor &c1, const Cursor &c2) {
                return OverlapSegment::compare(*c1, *c2);
            }) - cursors.begin();
            cursors[minIndex]++;
            
            // Stop the whole loop when any one iterator reaches the end
            if (cursors[minIndex] == candidates[minIndex].end())
                break;
        }
    }
}

template <class Hasher>
bool OverlapFinder<Hasher>::checkOverlap(vector<Residue *>::iterator begin, vector<Residue *>::iterator end, const OverlapSegment &candidate) const {
    if (!_verifier)
        return true;
                
    // Build lists of residues for each segment
    vector<Residue *> segment1(begin, end);
    
    vector<Residue *> chainResidues = candidate.structure->getChainByID(candidate.chainID)->getResidues();
    vector<Residue *> segment2(chainResidues.begin() + candidate.offset, chainResidues.begin() + candidate.offset + _segmentLength);
    
    // Verify
    return _verifier->verify(segment1, segment2);
}
