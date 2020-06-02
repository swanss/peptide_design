#include "clusterutils.h"
#include "Util.h"
#include "Fragments.h"
#include <unordered_set>

// FragmentFetcher

bool FragmentFetcher::hasNext() {
     if (!nextResult) {
        makeNextResult();
     }
     return nextResult != nullptr;
}

pair<AtomPointerVector, FragmentInfo *> FragmentFetcher::next() {
    if (!nextResult) {
        makeNextResult();
    }
    auto valueToReturn = make_pair(nextResult->first, nextResult->second->copy());
    destroyNextResult();
    return valueToReturn;
}

void FragmentFetcher::skip() {
    if (!nextResult) {
        makeNextResult();
    }
    destroyNextResult();
}

void FragmentFetcher::destroyNextResult() {
    // We own the memory for the pair and the fragment, so destroy both
    if (nextResult == nullptr) {
        return;
    }
    if (nextResult->second != nullptr) {
        delete nextResult->second;
    }
    delete nextResult; 
    nextResult = nullptr;
}

void formatStructure(Structure *s) {
    removeSideChains(*s);
    s->renumber(0, 0);
}

Structure *getFormattedStructure(StructureCache &cache, string name) {
    Structure *s;
    if (!cache.hasStructure(name)) {
        s = cache.getStructure(name);
        bool isValid = true;
        for (Residue *res: s->getResidues()) {
            if (!RotamerLibrary::hasFullBackbone(res)) {
                isValid = false;
                break;
            }
        }
        if (!isValid) {
            cache.removeStructure(name);
            return nullptr;
        } 

        formatStructure(s);
    } else {
        s = cache.getStructure(name);
    }
    return s;
}

// PairFragmentFetcher

PairFragmentFetcher::PairFragmentFetcher(string pdbsPath, string contactsPath, int numFlank, int splitCount): _numFlank(numFlank), _structureCache(), _splitCount(splitCount) {
    _structureCache.storeAtoms = true;
    pdbs = MstUtils::fileToArray(pdbsPath);
    contactFiles = MstUtils::fileToArray(contactsPath);
    _splitSize = pdbs.size() / _splitCount;
}

PairFragmentFetcher::~PairFragmentFetcher() {
    if (currentStructure != nullptr) {
        delete currentStructure;
        currentStructure = nullptr;
    }
}

void PairFragmentFetcher::makeNextResult() {
    if (nextResult != nullptr || _splitEnded)
        return;

    if (currentAPVIdx >= currentAPVs.size()) {
        // Time to get APVs from another structure
        currentAPVs.clear();
        currentAPVIdx = 0;

        while (currentAPVs.empty() && currentStructureIdx < (int)pdbs.size() - 1) {
            if (_splitCount != 1 && currentStructureIdx >= 0 && currentStructureIdx + 1 % _splitSize == 0) {
                currentStructureIdx++;
                cout << "End of mini-batch" << endl;
                _splitEnded = true;
                break;
            } else if (_splitCount == 1 || currentStructureIdx % _splitSize != 0) {
                currentStructureIdx++;
            }
            cout << "Loading a new structure, " << currentStructureIdx << endl;
            loadAPVsFromStructure(currentStructureIdx);
        } 
    }

    if (currentAPVs.empty()) {
        // Reached the end of the pdb list with no more results
        destroyNextResult();
        return;
    }

    // Fill nextResult with appropriate APV and info
    AtomPointerVector apv = currentAPVs[currentAPVIdx];
    if (apv.size() == numAtoms()) {
        nextResult = new pair<AtomPointerVector, FragmentInfo *>(apv, new PairFragmentInfo(currentStructureIdx, apv[0]->getIndex(), apv[numAtoms() / 2]->getIndex()));
        currentAPVIdx++;
    } else {
        cout << "Expected APV to have " << numAtoms() << " atoms, got " << apv.size() << endl;
        currentAPVIdx++;
        makeNextResult();
    }
}

bool PairFragmentFetcher::loadAPVsFromStructure(int structureIdx) {
    if (!MstSys::fileExists(pdbs[structureIdx]) || !MstSys::fileExists(contactFiles[structureIdx])) return false;
    
    Structure *structure = getFormattedStructure(_structureCache, pdbs[structureIdx]);
    if (!structure) return false;
    vector<Residue*> residues = structure->getResidues();
    contactList cl;
    
    // Keep a copy of the structure in case the structure cache deletes this one
    if (currentStructure != nullptr) {
        delete currentStructure;
    }
    currentStructure = new Structure(*structure);

    vector<string> lines = MstUtils::fileToArray(contactFiles[structureIdx]);
    for (int j = 0; j < lines.size(); j++) {
        vector<string> tmp = MstUtils::split(lines[j]);
        int resiA = stoi(tmp[1]); // ??
        int resiB = stoi(tmp[2]); // ??
        cl.addContact(residues[resiA], residues[resiB], 1.0);
    }
    
    Fragmenter fragmenter(*currentStructure);
    // make these a parameter
    FragmentParams fragParams(_numFlank, true, false);
    fragmenter.setParams(fragParams);
    fragmenter.fragment(cl);
    set<Fragment> frags = fragmenter.getFragments();
    for (auto frag: frags) {
        if (frag.numSegments() == 2) { // may want to change
            vector<Residue*> fragResidues = frag.getResidues();
            AtomPointerVector apv = residuesToAtoms(fragResidues);
            if (apv.size() != numAtoms()) {
                continue;
            }
            currentAPVs.push_back(apv);

            vector<vector<Residue*>> segs = frag.getSegments();
            fragResidues = segs[1];
            fragResidues.insert(fragResidues.end(), segs[0].begin(), segs[0].end());
            AtomPointerVector revAPV = residuesToAtoms(fragResidues);
            currentAPVs.push_back(revAPV);
        }
    }

    return true;
}

void PairFragmentFetcher::reset() {
    currentStructureIdx = -1;
    currentAPVIdx = 0;
    currentAPVs.clear();
}

AtomPointerVector PairFragmentFetcher::getAPV(FragmentInfo *info) {
    PairFragmentInfo *pairInfo = dynamic_cast<PairFragmentInfo *>(info);

    Structure *s = getFormattedStructure(_structureCache, pdbs[pairInfo->structureIdx]);
    int atomIdx1 = pairInfo->atomIdx1;
    int atomIdx2 = pairInfo->atomIdx2;
    
    // The only reliable pieces of information we have are numFlank, the start
    // atom indexes of each fragment, and the atom and residue indexes in the 
    // structure. Therefore, find the residue index of the first atom in each
    // chain, then walk up numFlank * 2 + 1 residues.
    vector<Atom*> &atoms = _structureCache.getAtoms(pdbs[pairInfo->structureIdx]);
    AtomPointerVector apv(numAtoms());

    Residue *res = atoms[atomIdx1]->getResidue();
    int atomIndex = 0;
    vector<Residue *> allResidues = s->getResidues();
    for (int residueOffset = 0; residueOffset < _numFlank * 2 + 1; residueOffset++) {
        for (Atom *atom: allResidues[res->getResidueIndex() + residueOffset]->getAtoms()) {
            apv[atomIndex++] = atom;
        } 
    }
    res = atoms[atomIdx2]->getResidue();
    for (int residueOffset = 0; residueOffset < _numFlank * 2 + 1; residueOffset++) {
        for (Atom *atom: allResidues[res->getResidueIndex() + residueOffset]->getAtoms()) {
            apv[atomIndex++] = atom;
        } 
    }

    return apv;
}

Structure *PairFragmentFetcher::getFullStructure(FragmentInfo *info) {
    PairFragmentInfo *pairInfo = dynamic_cast<PairFragmentInfo *>(info);
    return getFormattedStructure(_structureCache, pdbs[pairInfo->structureIdx]);
}

Structure PairFragmentFetcher::getResultStructure(FragmentInfo *info) {
    PairFragmentInfo *pairInfo = dynamic_cast<PairFragmentInfo *>(info);

    Structure *s = getFormattedStructure(_structureCache, pdbs[pairInfo->structureIdx]);
    int atomIdx1 = pairInfo->atomIdx1;
    int atomIdx2 = pairInfo->atomIdx2;
    
    // The only reliable pieces of information we have are numFlank, the start
    // atom indexes of each fragment, and the atom and residue indexes in the 
    // structure. Therefore, find the residue index of the first atom in each
    // chain, then walk up numFlank * 2 + 1 residues.
    vector<Atom*> atoms = s->getAtoms();
    vector<Residue *> ret;

    Residue *res = atoms[atomIdx1]->getResidue();
    vector<Residue *> allResidues = s->getResidues();
    for (int residueOffset = 0; residueOffset < _numFlank * 2 + 1; residueOffset++) {
        ret.push_back(allResidues[res->getResidueIndex() + residueOffset]);
    }
    res = atoms[atomIdx2]->getResidue();
    for (int residueOffset = 0; residueOffset < _numFlank * 2 + 1; residueOffset++) {
        ret.push_back(allResidues[res->getResidueIndex() + residueOffset]); 
    }

    return residuesToStructure(ret);
}

int PairFragmentFetcher::numAtoms() {
    return 8 * (_numFlank * 2 + 1);
}

FragmentInfo *PairFragmentFetcher::makeFragmentInfo(string infoString) {
    return new PairFragmentInfo(infoString);
}

// PairFragmentFetcher

SingleFragmentFetcher::SingleFragmentFetcher(StructuresBinaryFile *binaryFile, int numResidues, string chainID, int splitCount): _binaryFile(binaryFile), _cacheBinaryFile(new StructuresBinaryFile(*binaryFile)), _numResidues(numResidues), chainID(chainID), _splitCount(splitCount) {
    _structureCache = new StructureCache(_cacheBinaryFile, 10000000);
    _structureCache->storeAtoms = true;
    _binaryFile->scanFilePositions();
    _binaryFile->insertStructureNames(pdbs);
    _splitSize = pdbs.size() / _splitCount;
}

SingleFragmentFetcher::SingleFragmentFetcher(StructuresBinaryFile *binaryFile, StructureCache *cache, int numResidues, string chainID, int splitCount): _binaryFile(binaryFile), _cacheBinaryFile(new StructuresBinaryFile(*binaryFile)), _numResidues(numResidues), _structureCache(cache), chainID(chainID), _splitCount(splitCount), ownsStructureCache(false) {
    _structureCache->storeAtoms = true;
    _binaryFile->scanFilePositions();
    _binaryFile->insertStructureNames(pdbs);
    _splitSize = pdbs.size() / _splitCount;
}

SingleFragmentFetcher::~SingleFragmentFetcher() {
    if (currentStructure != nullptr) {
        delete currentStructure;
        currentStructure = nullptr;
    }
    if (_cacheBinaryFile != nullptr)
        delete _cacheBinaryFile;
    if (ownsStructureCache && _structureCache != nullptr) {
        delete _structureCache;
        _structureCache = nullptr;
    }
}

void SingleFragmentFetcher::makeNextResult() {
    if (nextResult != nullptr || _splitEnded)
        return;

    if (currentAPVIdx >= currentAPVs.size()) {
        // Time to get APVs from another structure
        currentAPVs.clear();
        currentAPVIdx = 0;

        while (currentAPVs.empty() && currentStructureIdx < (int)pdbs.size() - 1) {
            if (_splitCount != 1 && currentStructureIdx >= 0 && currentStructureIdx + 1 % _splitSize == 0) {
                currentStructureIdx++;
                cout << "End of mini-batch" << endl;
                _splitEnded = true;
                break;
            } else if (_splitCount == 1 || currentStructureIdx % _splitSize != 0) {
                currentStructureIdx++;
            }
            if (pdbs.size() > 100000) {
                if (currentStructureIdx % 1000 == 0)
                    cout << "Loading a new structure, " << currentStructureIdx << endl;
            } else
                cout << "Loading a new structure, " << currentStructureIdx << endl;
            loadAPVsFromStructure(currentStructureIdx);
        } 
    }

    if (currentAPVs.empty()) {
        // Reached the end of the pdb list with no more results
        destroyNextResult();
        return;
    }

    // Fill nextResult with appropriate APV and info
    AtomPointerVector &apv = currentAPVs[currentAPVIdx];
    if (apv.size() == numAtoms()) {
        nextResult = new pair<AtomPointerVector, FragmentInfo *>(apv, new SingleFragmentInfo(currentStructureIdx, apv[0]->getIndex()));
        currentAPVIdx++;
    } else {
        cout << "Expected APV to have " << numAtoms() << " atoms, got " << apv.size() << endl;
        currentAPVIdx++;
        makeNextResult();
    }
}

bool SingleFragmentFetcher::loadAPVsFromStructure(int structureIdx) {
    Structure *structure;
    
    if (_binaryFile != nullptr) {
        if (!_binaryFile->hasNext())
            return false;
        structure = _binaryFile->next();
        formatStructure(structure);
        MstUtils::assert(pdbs[structureIdx] == structure->getName(), "Names don't match"); 
    } else {
        if (!MstSys::fileExists(pdbs[structureIdx])) return false;
        structure = getFormattedStructure(*_structureCache, pdbs[structureIdx]);
    }

    if (!structure) return false;
    
    // Keep a copy of the structure in case the structure cache deletes this one
    if (currentStructure != nullptr) {
        delete currentStructure;
    }
    currentStructure = new Structure(*structure);
    currentAPVs.clear();

    // Load all _numResidues-length segments from contiguous chains
    for (int chainIdx = 0; chainIdx < currentStructure->chainSize(); chainIdx++) {
        Chain &chain = currentStructure->getChain(chainIdx);
        if (!chainID.empty() && chainID != chain.getID())
            continue;
        
        for (int resBase = 0; resBase < chain.residueSize() - _numResidues + 1; resBase++) {
            AtomPointerVector apv;
            for (int i = 0; i < _numResidues; i++) {
                vector<Atom *> atoms = chain.getResidue(resBase + i).getAtoms(); 
                for (Atom *a: atoms)
                    apv.push_back(a);
            }
            if (apv.size() != numAtoms())
                continue;
            currentAPVs.push_back(apv);
        }
    }

    return true;
}

void SingleFragmentFetcher::reset() {
    if (currentStructure != nullptr)
        delete currentStructure;
    currentStructure = nullptr;
    destroyNextResult();
    currentStructureIdx = -1;
    if (_binaryFile != nullptr)
        _binaryFile->reset();
    //_structureCache.clear();
}

AtomPointerVector SingleFragmentFetcher::getAPV(FragmentInfo *info) {
    SingleFragmentInfo *fragInfo = dynamic_cast<SingleFragmentInfo *>(info);
    int atomIdx = fragInfo->atomIdx;

    MstUtils::assert(fragInfo->structureIdx < pdbs.size(), "Requested structure index out of bounds");
    Structure *s = getFormattedStructure(*_structureCache, pdbs[fragInfo->structureIdx]);
    vector<Atom *> &atoms  = (*_structureCache).getAtoms(pdbs[fragInfo->structureIdx]);
    AtomPointerVector apv(numAtoms());

    Residue *res = atoms[atomIdx]->getResidue();
    int atomIndex = 0;
    vector<Residue *> allResidues = s->getResidues();
    for (int residueOffset = 0; residueOffset < _numResidues; residueOffset++) {
        for (Atom *atom: allResidues[res->getResidueIndex() + residueOffset]->getAtoms()) {
            apv[atomIndex++] = atom;
        } 
    }
    return apv;
}

Structure *SingleFragmentFetcher::getFullStructure(FragmentInfo *info) {
    SingleFragmentInfo *fragInfo = dynamic_cast<SingleFragmentInfo *>(info);
    return getFormattedStructure(*_structureCache, pdbs[fragInfo->structureIdx]);
}

Structure SingleFragmentFetcher::getResultStructure(FragmentInfo *info) {
    SingleFragmentInfo *fragInfo = dynamic_cast<SingleFragmentInfo *>(info);

    Structure *s = getFormattedStructure(*_structureCache, pdbs[fragInfo->structureIdx]);
    int atomIdx = fragInfo->atomIdx;
    
    vector<Atom*> atoms = s->getAtoms();
    vector<Residue *> ret;

    Residue *res = atoms[atomIdx]->getResidue();
    vector<Residue *> allResidues = s->getResidues();
    for (int residueOffset = 0; residueOffset < _numResidues; residueOffset++) {
        ret.push_back(allResidues[res->getResidueIndex() + residueOffset]);
    }
 
    return residuesToStructure(ret);
}

int SingleFragmentFetcher::numAtoms() {
    return 4 * _numResidues;
}

FragmentInfo *SingleFragmentFetcher::makeFragmentInfo(string infoString) {
    return new SingleFragmentInfo(infoString);
}

// BatchWorkerFragmentFetcher

void BatchWorkerFragmentFetcher::makeNextResult() {
    while (_fragmentIndex % _numWorkers != _workerIndex && _fetcher->hasNext()) {
        _fetcher->skip();
        _fragmentIndex++;
    }
    _fragmentIndex++;
    if (_fetcher->hasNext())
        nextResult = new pair<AtomPointerVector, FragmentInfo *>(_fetcher->next());
    else
        nextResult = nullptr;
}

void BatchWorkerFragmentFetcher::reset() {
    _fetcher->reset();
    _fragmentIndex = 0;
}

AtomPointerVector BatchWorkerFragmentFetcher::getAPV(FragmentInfo *info) {
    return _fetcher->getAPV(info);
}

int BatchWorkerFragmentFetcher::numAtoms() {
    return _fetcher->numAtoms();
}

FragmentInfo *BatchWorkerFragmentFetcher::makeFragmentInfo(string infoString) {
    return _fetcher->makeFragmentInfo(infoString);
}

Structure *BatchWorkerFragmentFetcher::getFullStructure(FragmentInfo *info) {
    return _fetcher->getFullStructure(info); 
}

Structure BatchWorkerFragmentFetcher::getResultStructure(FragmentInfo *info) {
    return _fetcher->getResultStructure(info);
}

