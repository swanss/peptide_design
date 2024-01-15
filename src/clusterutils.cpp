#include "clusterutils.h"
#include "utilities.h"
#include "fragments.h"
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
        MstUtils::assertCond(pdbs[structureIdx] == structure->getName(), "Names don't match"); 
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

    MstUtils::assertCond(fragInfo->structureIdx < pdbs.size(), "Requested structure index out of bounds");
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

/*
 Loads each seed and finds the windows
 */
GreedyClusterer::GreedyClusterer(string seedBin_path, int _window_length) {
    window_length = _window_length;
    
    cout << "Loading seeds and finding windows..." << endl;
    StructureIterator* structIt = new StructureIterator(seedBin_path);
    
    while (structIt->hasNext()) {
        vector<Structure*> seedStructures = structIt->next();
        for (Structure* seed : seedStructures) {
            int num_windows = seed->residueSize() - window_length + 1;
            for (int i = 0; i < num_windows; i++) {
                seedWindowInfo info(seed->getName(),i);
                int unique_ID = seedWindows.size();
                seedWindows[unique_ID] = info;
                seedWindowsRev[info] = unique_ID;
                remainingSeedWindows.insert(unique_ID);
            }
        }
    }
    delete structIt;
    
    cout << "Constructing seed cache..." << endl;
    seedBin = new StructuresBinaryFile(seedBin_path);
    seedStructures = StructureCache(seedBin,10000);
    
    // initialize overlaps variable
    for (int ID : remainingSeedWindows) {
        seedWindowOverlaps[ID] = set<int>();
    }
    numTotalSeedWindows = remainingSeedWindows.size();
    cout << "Found " << seedBin->structureCount() << " seeds and "<< seedWindows.size() << " windows" << endl;
};

void GreedyClusterer::addOverlapInfo(FuseCandidateFile file) {
    vector<FuseCandidate> candidateBatch;
    int count = 0;
    while ((candidateBatch = file.read(1)).size() > 0) {
        FuseCandidate cand = candidateBatch[0];
        if (cand.overlapSize != window_length) MstUtils::error("Overlap size in fuse candidate file does not match the window size that was specified during construction");
        seedWindowInfo window_A(cand.file1,cand.overlapPosition1);
        seedWindowInfo window_B(cand.file2,cand.overlapPosition2);
        
        // verify that the overlapping windows are in the set we're trying to cover
        if (seedWindowsRev.count(window_A) == 0) MstUtils::error("Window provided in overlap not in the set that will be covered: "+window_A.getName(),"GreedyClusterer::addOverlapInfo");
        if (seedWindowsRev.count(window_B) == 0) MstUtils::error("Window provided in overlap not in the set that will be covered: "+window_B.getName(),"GreedyClusterer::addOverlapInfo");
        
        // store the overlap information
        seedWindowOverlaps[seedWindowsRev[window_A]].insert(seedWindowsRev[window_B]);
        seedWindowOverlaps[seedWindowsRev[window_B]].insert(seedWindowsRev[window_A]);
        count++;
    }
    cout << "Added " << count << " overlaps." << endl;
};

void GreedyClusterer::performClustering(mstreal max_coverage, int max_clusters) {
    int numCoveredReq = max_coverage*numTotalSeedWindows;
    int numCovered = 0;
    cout << "Will terminate when " << numCoveredReq << " are covered";
    if (max_clusters > 0) cout << " or when " << max_clusters << " clusters identified" << endl;
    else cout << endl;

    int nCluster = 0;
    while ((numCovered < numCoveredReq)||((max_clusters > 0)&&(nCluster <= max_clusters))) {
        
        // iterate over the clusters to find the largest
        int max_cluster_ID = *remainingSeedWindows.begin();
        for (int ID : remainingSeedWindows) {
            if (seedWindowOverlaps[ID].size() > seedWindowOverlaps[max_cluster_ID].size()) max_cluster_ID = ID;
        }
        
        // Store the cluster center as well as the members
        clusters.push_back(pair<int,set<int>>(max_cluster_ID,seedWindowOverlaps[max_cluster_ID]));
        
        /*
        In order to prepare for the next iteration ALL information pertaining to members of this cluster
        must be removed. This includes the following:
        A) All covered windows must be removed from remainingSeedWindows
        B) All overlaps to the covered windows from seedWindowOverlaps
         
        Since the overlap information is duplicated, we will need to delete each overlap twice. Consider
        the following overlap data:
         
          A  B  C
        A X  X
        B X  X  X
        C    X  X
         
        If the cluster around seed A is removed, we will need to delete row A since it contains all of the overlaps to
        A, row B since it overlaps A and is therefore in the cluster, and entry B from row C since B is covered and 
        that overlap is no longer considered.
         
          A  B  C
        A -  -
        B -  -  -
        C    -  X
        
        */
        
        // remove the cluster center from the set of remaining windows
        if (remainingSeedWindows.count(max_cluster_ID) == 0) MstUtils::error("Cannot delete "+MstUtils::toString(max_cluster_ID)+" not in the remainingSeedWindows","GreedyClusterer::performClustering");
        remainingSeedWindows.erase(max_cluster_ID);
        for (int ID_b : seedWindowOverlaps[max_cluster_ID]) {
            // remove cluster member from the set of remaining windows
            if (remainingSeedWindows.count(ID_b) == 0) MstUtils::error("Cannot delete "+MstUtils::toString(ID_b)+" not in the remainingSeedWindows","GreedyClusterer::performClustering");
            remainingSeedWindows.erase(ID_b);
            
            for (int ID_c : seedWindowOverlaps[ID_b]) {
                // remove the overlap between cluster member and other window
                if (seedWindowOverlaps[ID_c].count(ID_b) == 0) MstUtils::error("Cannot delete overlap "+MstUtils::toString(ID_b)+" not in the seedWindowOverlaps","GreedyClusterer::performClustering");
                seedWindowOverlaps[ID_c].erase(ID_b);
            }
            // remove all overlaps to the cluster member
            if (seedWindowOverlaps.count(ID_b) == 0) MstUtils::error("Cannot delete all overlaps to "+MstUtils::toString(ID_b)+" not in the seedWindowOverlaps","GreedyClusterer::performClustering");
            seedWindowOverlaps.erase(ID_b);
        }
        // remove all overlaps to the cluster center
        if (seedWindowOverlaps.count(max_cluster_ID) == 0) MstUtils::error("Cannot delete "+MstUtils::toString(max_cluster_ID)+" not in the seedWindowOverlaps","GreedyClusterer::performClustering");
        seedWindowOverlaps.erase(max_cluster_ID);
        
        // check coverage progress
        numCovered = numTotalSeedWindows - remainingSeedWindows.size();
        
        cout << "Cluster: " << nCluster << " with " << clusters.back().second.size()+1 << " members. Overall " << numCovered << " of the seeds are covered" << endl;
        nCluster++;
    }
    cout << "Done clustering. In total there are " << clusters.size() << " clusters" << endl;
}

void GreedyClusterer::writeClusterInfo(string outDir, Chain* peptide, bool pdbs, bool verbose) {
    map<int,pair<int,mstreal>> cluster2peptideWindow; // Maps the cluster ID to the N-terminal residue of the overlapping peptide window
    if (peptide != nullptr) {
        // If a peptide chain is provided, try to map the clusters to the peptide chain
        RMSDCalculator rmsdCalc;
        
        vector<Atom*> peptideAtoms;
        for (Residue* R : peptide->getResidues()) {
            if (RotamerLibrary::hasFullBackbone(R)) {
                vector<Atom*> peptideResAtoms = RotamerLibrary::getBackbone(R);
                peptideAtoms.insert(peptideAtoms.end(),peptideResAtoms.begin(),peptideResAtoms.end());
            } else MstUtils::error("Residue "+R->getChainID()+MstUtils::toString(R->getNum())+" missing backbone atoms");
        }
        cout << "Peptide chain with " << peptide->residueSize() << " residues and " << peptideAtoms.size() << " atoms" << endl;
        for (int i = 0; i < clusters.size(); i++) {
            Structure cluster_rep_seed = getSeedWindowStructure(seedWindows[clusters[i].first]);
            vector<Atom*> cluster_rep_seed_atoms = cluster_rep_seed.getAtoms();
            for (int pos = 0; pos < peptide->residueSize() - window_length + 1; pos++) {
                vector<Atom*> peptide_window(peptideAtoms.begin()+(pos*4),peptideAtoms.begin()+(pos+window_length)*4);
                mstreal rmsd = rmsdCalc.rmsd(peptide_window,cluster_rep_seed_atoms);
                if (pos == 0) cluster2peptideWindow[i] = pair<int,mstreal>(pos,rmsd);
                else if (rmsd < cluster2peptideWindow[i].second) cluster2peptideWindow[i] = pair<int,mstreal>(pos,rmsd);
            }
        }
    }
    
    fstream out;
    string filename = outDir + "/cluster_info.csv";
    MstUtils::openFile(out,filename,fstream::out);
    out << "cluster_number,cluster_size,cluster_representative,overlapping_peptide_window,rmsd" << endl;
    for (int i = 0; i < clusters.size(); i++) {
        string seed_window_name = seedWindows[clusters[i].first].getName();
        out << i << ",";
        out << clusters[i].second.size()+1 << ",";
        out << seed_window_name << ",";
        out << ((peptide != nullptr) ? MstUtils::toString(cluster2peptideWindow[i].first) : "") << ",";
        out << ((peptide != nullptr) ? MstUtils::toString(cluster2peptideWindow[i].second) : "") << endl;
    }
    out.close();
    
    fstream all_seed_windows;
    filename = outDir + "/cluster_members.csv";
    if (verbose) MstUtils::openFile(all_seed_windows,filename,fstream::out);
    
    if (pdbs) {
        for (int i = 0; i < clusters.size(); i++) {
            string cluster_centroid_name = seedWindows[clusters[i].first].getName();
            Structure cluster_centroid = getSeedWindowStructure(seedWindows[clusters[i].first]);
            string pdb_name = outDir+MstUtils::toString(i)+"_"+MstUtils::toString(0)+"_"+cluster_centroid_name+".pdb";
            cluster_centroid.writePDB(pdb_name);
            if (verbose) all_seed_windows << cluster_centroid_name;
            int j = 1;
            for (auto it = clusters[i].second.begin(); it != clusters[i].second.end(); it++) {
                if (*it == clusters[i].first) continue; //duplicate of the representative seed
                string seed_window_name = seedWindows[*it].getName();
                Structure seed_window = getSeedWindowStructure(seedWindows[*it]);
                string pdb_name = outDir+MstUtils::toString(i)+"_"+MstUtils::toString(j)+"_"+seed_window_name+".pdb";
                seed_window.writePDB(pdb_name);
                if (verbose) all_seed_windows << "," << seed_window_name;
                j++;
            }
            if (verbose) all_seed_windows << endl;
        }
    }
    all_seed_windows.close();
}
