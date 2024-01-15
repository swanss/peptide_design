#include <stdio.h>
#include <stack>
#include "seedgraph.h"
#include "findfuseable.h"
#include "utilities.h"

// SeedGraph

SeedGraph::SeedGraph(bool adjSameResidues, StructureCache *sCache, bool centerOnly): adjSameResidues(adjSameResidues), centerOnly(centerOnly) {
    if (adjSameResidues == true and centerOnly == true) MstUtils::error("centerOnly cannot be applied while adjSameResidues is true");
    if (sCache == nullptr)
        structures = new StructureCache();
    else {
        structures = sCache;
        ownsStructureCache = false;
    }
}

SeedGraph::SeedGraph(string adjacencyPath, bool adjSameResidues, StructureCache *sCache, string pdbPrefix, bool centerOnly): adjSameResidues(adjSameResidues), centerOnly(centerOnly) {
    if (adjSameResidues == true and centerOnly == true) MstUtils::error("centerOnly cannot be applied while adjSameResidues is true");
    if (sCache == nullptr)
        structures = new StructureCache(pdbPrefix);
    else {
        structures = sCache;
        ownsStructureCache = false;
    }
    read(adjacencyPath, pdbPrefix);
}

SeedGraph::~SeedGraph() {
    if (structures != nullptr && ownsStructureCache) {
        delete structures;
        structures = nullptr;
    }
}

void SeedGraph::load(FuseCandidateFile neighborsFile, string pdbPrefix, Structure *nearStructure, float inset, FuseCandidateFile *usedCandidatesFile) {
    if (structures->getPDBPrefix().size() == 0)
        structures->setPDBPrefix(pdbPrefix);
    vector<FuseCandidate> candidateBatch;
    FuseCandidateFinder fuseHelper;
    int fuseIndex = 0;
    while ((candidateBatch = neighborsFile.read(1)).size() > 0) {
        FuseCandidate cand = candidateBatch[0];
        
        Structure *structure1 = structures->getStructure(cand.file1, pdbPrefix);
        Structure *structure2 = structures->getStructure(cand.file2, pdbPrefix);
        
        // Check if either structure overlaps with the bounding box of nearStructure, if applicable
        if (nearStructure != nullptr) {
            vector<Structure *> toCheck = { structure1, structure2 };
            if (fuseHelper.findNearbySeeds(nearStructure, toCheck, inset).size() < 2)
                continue;
            if (usedCandidatesFile != nullptr) {
                // Write out that we used this candidate
                usedCandidatesFile->write(cand);
            }
        }
        
        Chain *chain1 = structure1->getChainByID(cand.chain1);
        Chain *chain2 = structure2->getChainByID(cand.chain2);
        
        addAdjacencies(chain1, chain2, cand.overlapPosition1, cand.overlapPosition2, cand.overlapSize);
        addAdjacencies(chain2, chain1, cand.overlapPosition2, cand.overlapPosition1, cand.overlapSize);
    }
}

void SeedGraph::load(SeedListFile structuresFile, string pdbPrefix, Structure *nearStructure, float inset) {
    if (structures->getPDBPrefix().size() == 0)
        structures->setPDBPrefix(pdbPrefix);
    FuseCandidateFinder fuseHelper;
    vector<string> filePaths = structuresFile.read().first;
    for (string path: filePaths) {
        Structure *structure = structures->getStructure(path, pdbPrefix);
        
        // Check if either structure overlaps with the bounding box of nearStructure, if applicable
        if (nearStructure != nullptr) {
            vector<Structure *> toCheck = { structure };
            if (fuseHelper.findNearbySeeds(nearStructure, toCheck, inset).size() < 1)
                continue;
        }
        
        Chain *chain = structure->getChainByID(seed_chain_ID);
        
        addAdjacencies(chain, chain, 0, 0, 0);
    }
}

void SeedGraph::load(StructuresBinaryFile *binaryFile, Structure *nearStructure, float inset) {
    FuseCandidateFinder fuseHelper;
    binaryFile->reset();
    while (binaryFile->hasNext()) {
        Structure *structure = binaryFile->next();
        
        // Check if either structure overlaps with the bounding box of nearStructure, if applicable
        if (nearStructure != nullptr) {
            vector<Structure *> toCheck = { structure };
            if (fuseHelper.findNearbySeeds(nearStructure, toCheck, inset).size() < 1)
                continue;
        }
        
        Chain *chain = structure->getChainByID(seed_chain_ID);
        
        addAdjacencies(chain, chain, 0, 0, 0);
    }
}

void SeedGraph::loadCache() {
    //check if the cache has preloaded all of the seeds in the binary file
    if (!structures->isPreloaded()) structures->preloadFromBinaryFile();
    if (!structures->belowCapacity()) MstUtils::error("Unable to load from cache as the capacity is lower than the number of seeds in the binary file","SeedGraph::loadCache()");
    for (auto it = structures->begin(); it != structures->end(); ++it) {
        Structure *structure = *it;
        Chain *chain = structure->getChainByID(seed_chain_ID);
        
        addAdjacencies(chain, chain, 0, 0, 0);
    }
}


void SeedGraph::load(ClusterTree &clusterTree, int numResOverlap, mstreal rmsdCutoff, mstreal minCosAngle) {
    FragmentFetcher *fetcher = clusterTree.getFragmentFetcher();
    FuseCandidateFinder overlapDetector(numResOverlap, FuseCandidateSearchMode::general, rmsdCutoff);
    overlapDetector.minCosAngle = minCosAngle;

    ClusterPairIterator it(clusterTree);
    unsigned long long numOverlapsFound = 0;
    unsigned long long numNodesSkipped = 0;
    unsigned long long numNodesAutoAdded = 0;
    unsigned long long pairIdx = 0;

    ClusterNode<FragmentInfo> *currentFirstNode = nullptr;
    long numSecondNodesTested = 0;
    long totalSecondNodesTested = 0;
    long totalFirstNodes = 0;

    while (it.hasNext()) {
        if (++pairIdx % 1000000 == 0)
            cout << "Pair " << pairIdx << " - " << numOverlapsFound << " overlaps, " << numNodesSkipped << " skipped, " << numNodesAutoAdded << " auto-added" << endl;
        auto pairItem = it.next();

        // Progress tracking
        if (currentFirstNode != pairItem.first) {
            if (currentFirstNode) {
                totalSecondNodesTested += numSecondNodesTested;
                totalFirstNodes++;
                if (totalFirstNodes % 100 == 0) {
                    cout << totalFirstNodes << " first nodes tested, avg. " << totalSecondNodesTested / 100.0 << " comparisons each" << endl;
                    totalSecondNodesTested = 0;
                }
            }
            currentFirstNode = pairItem.first;
            numSecondNodesTested = 0;
        }

        if (pairItem.first == pairItem.second) {
            if (pairItem.first->subtreeRadius <= rmsdCutoff) {
                // All children of this node will be within the RMSD cutoff!
                it.skipSecondSubtree();

                AtomPointerVector apv1 = fetcher->getAPV(pairItem.first->getItem());
                Residue *res1 = apv1[0]->getParent();
                Chain *chain1 = res1->getParent();
                int pos1 = res1->getResidueIndexInChain();

                ClusterIterator childIt(pairItem.second);
                int numAdded = 0;
                while (childIt.hasNext()) {
                    auto node = childIt.next();
                    if (node == pairItem.first)
                        continue;
                    AtomPointerVector apv3 = fetcher->getAPV(node->getItem());
                    Residue *res3 = apv3[0]->getParent();
                    Chain *chain3 = res3->getParent();
                    int pos3 = res3->getResidueIndexInChain();
                    addAdjacencies(chain1, chain3, pos1, pos3, numResOverlap);
                    addAdjacencies(chain3, chain1, pos3, pos1, numResOverlap);
                    numOverlapsFound++;
                    numAdded++;
                }
                numNodesAutoAdded++;
            }
            continue;
        }

        AtomPointerVector apv1 = fetcher->getAPV(pairItem.first->getItem());
        AtomPointerVector apv2 = fetcher->getAPV(pairItem.second->getItem());

        // In order to completely eliminate the second subtree, we need to lower
        // bound the overlap between these two nodes relative to the radius of
        // the second subtree.
        mstreal rmsd = overlapDetector.quickRMSD(apv1, apv2, rmsdCutoff + pairItem.second->subtreeRadius);
        if (rmsd < 1000) {
            numSecondNodesTested++;
            rmsd = overlapDetector.atomsFormOverlap(apv1, apv2, rmsdCutoff);
            Residue *res1 = apv1[0]->getParent();
            Chain *chain1 = res1->getParent();
            int pos1 = res1->getResidueIndexInChain();

            if (rmsd <= rmsdCutoff) {
                // Add overlaps
                Residue *res2 = apv2[0]->getParent();
                Chain *chain2 = res2->getParent();
                int pos2 = res2->getResidueIndexInChain();
                addAdjacencies(chain1, chain2, pos1, pos2, numResOverlap);
                addAdjacencies(chain2, chain1, pos2, pos1, numResOverlap);
                numOverlapsFound++;
            }
            if (rmsd + pairItem.second->subtreeRadius <= rmsdCutoff) {
                // All children of this node will have an overlap!
                it.skipSecondSubtree();

                ClusterIterator childIt(pairItem.second);
                int numAdded = 0;
                while (childIt.hasNext()) {
                    auto node = childIt.next();
                    if (node == pairItem.second)
                        continue;
                    AtomPointerVector apv3 = fetcher->getAPV(node->getItem());
                    Residue *res3 = apv3[0]->getParent();
                    Chain *chain3 = res3->getParent();
                    int pos3 = res3->getResidueIndexInChain();
                    addAdjacencies(chain1, chain3, pos1, pos3, numResOverlap);
                    addAdjacencies(chain3, chain1, pos3, pos1, numResOverlap);
                    numOverlapsFound++;
                    numAdded++;
                }
                numNodesAutoAdded++;
            }
        } else {
            it.skipSecondSubtree();
            numNodesSkipped++;
        }
    }
    cout << numNodesSkipped << " nodes skipped, " << numNodesAutoAdded << " auto-added" << endl;
}

void SeedGraph::read(string adjacencyPath, string pdbPrefix, SeedGraph *convertToAdjSameResidue) {
    if (structures->getPDBPrefix().size() == 0)
        structures->setPDBPrefix(pdbPrefix);

    if (convertToAdjSameResidue != nullptr) {
        //SeedGraph tempGraph(adjacencyPath, false, structures, pdbPrefix);
        
        for (auto pair: convertToAdjSameResidue->adjacencies) {
            addResidue(pair.first);
            
            // Iterate over edges going forward, then backward
            for (Residue *next: pair.second) {
                for (Residue *prev: convertToAdjSameResidue->reverseAdjacencies[next]) {
                    if (prev != pair.first) {
                        addEdge(pair.first, prev);
                        addEdge(prev, pair.first);
                    }
                }
            }
            
            // Iterate over edges going backward, then forward
            for (Residue *prev: convertToAdjSameResidue->reverseAdjacencies[pair.first]) {
                for (Residue *next: convertToAdjSameResidue->adjacencies[prev]) {
                    if (next != pair.first) {
                        addEdge(pair.first, next);
                        addEdge(next, pair.first);
                    }
                }
            }
        }
        
        return;
    }
    
    ifstream readstream(adjacencyPath);
    if (!readstream.is_open()) {
        cerr << "couldn't open input stream (SeedGraph::read)" << endl;
        return;
    }
    
    string line;
    while (getline(readstream, line)) {
        vector<string> comps = splitString(line, ",");
        Residue *baseResidue = getResidueFromFile(comps[0]);
        addResidue(baseResidue);
        for (int i = 1; i < comps.size(); i++) {
            Residue *otherResidue = getResidueFromFile(comps[i]);
            addEdge(baseResidue, otherResidue);
        }
    }
}

Residue* SeedGraph::getResidueFromFile(string residueID, bool loadIfNeeded) {
    vector<string> comps = splitString(residueID, ":");
    if (!loadIfNeeded && !structures->hasStructure(comps[0]))
        return nullptr;
    
    Chain *chain = structures->getStructure(comps[0])->getChainByID(seed_chain_ID);
    return &chain->getResidue(atoi(comps[1].c_str()));
}

Structure* SeedGraph::getStructureFromFile(string structureName, bool loadIfNeeded) {
    if (!loadIfNeeded && !structures->hasStructure(structureName)) return nullptr;
    else {
        Structure *structure = structures->getStructure(structureName);
        return structure;
    }
}

void SeedGraph::print() {
    cout << "{" << endl;
    for (auto elem: adjacencies) {
        cout << "\t" << elem.first->getName() << elem.first->getResidueIndexInChain() << ": [ ";
        for (Residue *res: elem.second) {
            cout << res->getName() << res->getResidueIndexInChain() << " ";
        }
        cout << "]" << endl;
    }
    cout << "}" << endl;
}

void SeedGraph::addAdjacencies(Chain *chain1, Chain *chain2, int overlap1, int overlap2, int overlapSize) {
    if (centerOnly and (overlapSize % 2 != 0)) MstUtils::error("if centerOnly provided, overlapSize must be even");
    
    for (int i = 0; i < chain1->residueSize(); i++) {
        Residue &res = chain1->getResidue(i);
        addResidue(&res);
        
        //check if this position falls within an overlap region
        bool overlap_region = chain1 != chain2 && i >= overlap1 && i < overlap1 + overlapSize;
        
        //add applicable edges
        if (adjSameResidues) {
            if (overlap_region) {
                // Add an edge to the corresponding residue on the other chain
                int otherIdx = i - overlap1 + overlap2;
                if (otherIdx >= 0 && otherIdx < chain2->residueSize()) {
                    Residue &otherRes = chain2->getResidue(otherIdx);
                    addEdge(&res, &otherRes);
                }
            }
        } else {
            if (overlap_region) {
                if (!centerOnly) {
                    // Add an edge to the next residue on the other chain
                    int otherIdx = i - overlap1 + overlap2 + 1;
                    if (otherIdx >= 0 && otherIdx < chain2->residueSize()) {
                        Residue &otherRes = chain2->getResidue(otherIdx);
                        addEdge(&res, &otherRes);
                    }
                }
                else if (i == overlap1 + (overlapSize / 2) - 1) {
                    // edge will cross the center of the overlap
                    // Add an edge to the next residue on the other chain
                    int otherIdx = i - overlap1 + overlap2 + 1;
                    if (otherIdx >= 0 && otherIdx < chain2->residueSize()) {
                        Residue &otherRes = chain2->getResidue(otherIdx);
                        addEdge(&res, &otherRes);
                    }
                }
               
            }
            // Add an edge to the next residue on this chain
            if (i < chain1->residueSize() - 1) {
                Residue &otherRes = chain1->getResidue(i + 1);
                addEdge(&res, &otherRes);
            }
        }
    }
}

void SeedGraph::addResidue(Residue *res) {
    adjacencies[res];
    reverseAdjacencies[res];
}

void SeedGraph::addEdge(Residue *res1, Residue *res2) {
    numEdges++;
    adjacencies[res1].insert(res2);
    reverseAdjacencies[res2].insert(res1);
}

unordered_set<Residue *> SeedGraph::getResidues() {
    unordered_set<Residue *> residues;
    for (auto it: adjacencies) {
        residues.insert(it.first);
    }
    return residues;
}

vector<SeedGraph> SeedGraph::subgraphs() {
    unordered_set<Residue *> seenResidues;
    vector<SeedGraph> result;
    
    for (auto it: adjacencies) {
        if (seenResidues.find(it.first) != seenResidues.end())
            continue;
        
        // Get the neighborhood of this residue
        SeedGraph subgraph = neighborhood(it.first);
        result.push_back(subgraph);
        
        for (Residue *res: subgraph.getResidues()) {
            MstUtils::assertCond(seenResidues.find(res) == seenResidues.end(), "seen residues shouldn't already contain residue from subgraph");
            seenResidues.insert(res);
        }
    }
    
    return result;
}

unordered_set<Residue *> SeedGraph::bidirectionalNeighbors(Residue *res) {
    unordered_set<Residue *> result;
    auto adj = adjacencies[res];
    result.insert(adj.begin(), adj.end());
    auto revAdj = reverseAdjacencies[res];
    result.insert(revAdj.begin(), revAdj.end());
    return result;
}

unordered_set<Residue *> SeedGraph::forwardNeighbors(Residue *res) {
    MstUtils::assertCond(!adjSameResidues);
    return adjacencies[res];
}

unordered_set<Residue *> SeedGraph::backwardNeighbors(Residue *res) {
    MstUtils::assertCond(!adjSameResidues);
    return reverseAdjacencies[res];
}

SeedGraph SeedGraph::neighborhood(Residue *residue) {
    SeedGraph newGraph(adjSameResidues, structures);
    
    deque<Residue *> residuesToCheck;
    unordered_set<Residue *> seenResidues;
    residuesToCheck.push_back(residue);
    seenResidues.insert(residue);

    while (residuesToCheck.size() > 0) {
        Residue *res = residuesToCheck.front();
        newGraph.addResidue(res);
        residuesToCheck.pop_front();
        
        // Add edges from this node
        for (Residue *adjacent: adjacencies[res]) {
            if (seenResidues.count(adjacent) != 0)
                continue;
            newGraph.addEdge(res, adjacent);
            residuesToCheck.push_back(adjacent);
            seenResidues.insert(adjacent);
        }
        // Add edges into this node
        for (Residue *adjacent: reverseAdjacencies[res]) {
            if (seenResidues.count(adjacent) != 0)
                continue;
            newGraph.addEdge(adjacent, res);
            residuesToCheck.push_back(adjacent);
            seenResidues.insert(adjacent);
        }
    }
    
    return newGraph;
}

SeedGraph SeedGraph::removing(SeedGraph &subgraph) {
    SeedGraph newGraph = *this;
    for (Residue *res: subgraph.getResidues()) {
        newGraph.adjacencies.erase(res);
        newGraph.reverseAdjacencies.erase(res);
    }
    return newGraph;
}

SeedGraph SeedGraph::unionWith(SeedGraph &graph) {
    SeedGraph newGraph = *this;
    for (Residue *res: graph.getResidues()) {
        newGraph.addResidue(res);
        for (Residue *other: graph.adjacencies[res])
            newGraph.addEdge(res, other);
        for (Residue *other: graph.reverseAdjacencies[res])
            newGraph.addEdge(other, res);
    }
    return newGraph;
}

size_t SeedGraph::seedSize() {
    unordered_set<string> seenPaths;
    for (auto it: adjacencies) {
        string strName = it.first->getStructure()->getName();
        if (seenPaths.count(strName) == 0) {
            seenPaths.insert(strName);
        }
    }
    return seenPaths.size();
}

SeedGraph SeedGraph::withAdjSameResidue() {
    SeedGraph newGraph(true, structures);
    int i = 0;
    for (auto pair: adjacencies) {
        if (adjacencies.size() > 10000 && i % 10000 == 1)
            cout << "Residue " << i - 1 << " out of " << adjacencies.size() << endl;
        newGraph.addResidue(pair.first);
        
        // Iterate over edges going forward, then backward
        for (Residue *next: pair.second) {
            for (Residue *prev: reverseAdjacencies[next]) {
                if (prev != pair.first) {
                    newGraph.addEdge(pair.first, prev);
                    newGraph.addEdge(prev, pair.first);
                }
            }
        }

        // Iterate over edges going backward, then forward
        for (Residue *prev: reverseAdjacencies[pair.first]) {
            for (Residue *next: adjacencies[prev]) {
                if (next != pair.first) {
                    newGraph.addEdge(pair.first, next);
                    newGraph.addEdge(next, pair.first);
                }
            }
        }
        
        i++;
    }
    
    MstUtils::assertCond(newGraph.residueSize() == residueSize(), "unequal sizes");
    return newGraph;
}

string SeedGraph::writeCodeForResidue(Residue *res) {
    string structureName = res->getStructure()->getName();
    string writeName;
    if (structures->getPDBPrefix().length() > 0)
        writeName = MstSystemExtension::relativePath(structureName, structures->getPDBPrefix());
    else
        writeName = MstSystemExtension::fileName(structureName);
    return writeName + ":" + to_string(res->getResidueIndexInChain());
}

void SeedGraph::write(string path) {
    ofstream ss(path, ios::out);
    if (!ss.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return;
    }

    // Write the adjacencies directly
    for (auto it: adjacencies) {
        Residue *res = it.first;
        ss << writeCodeForResidue(res);
        if (it.second.size() > 0)
            ss << ",";
        int i = 0;
        for (Residue *otherRes: it.second) {
            ss << writeCodeForResidue(otherRes);
            if (++i != it.second.size())
                ss << ",";
        }
        ss << endl;
    }
}

vector<Residue *> SeedGraph::similarContactResidues(Residue *residue, int flankLength, float rmsdCutoff, vector<float> *rmsds) {
    vector<Residue *> results;
    
    MstUtils::assertCond(adjacencies.count(residue) != 0, "residue not found in graph");
    
    // Get the flanking region of this residue
    vector<Residue *> chainResidues = residue->getChain()->getResidues();
    int resIdx = residue->getResidueIndexInChain();
    int chainSize = residue->getChain()->residueSize();
    int realFlankLeft = resIdx >= flankLength ? flankLength : resIdx;
    int realFlankRight = resIdx < chainSize - flankLength ? flankLength : chainSize - resIdx - 1;
    vector<Atom *> mainFlank;
    for (int i = resIdx - realFlankLeft; i < resIdx + realFlankRight; i++) {
        vector<Atom *> resAtoms = chainResidues[i]->getAtoms();
        mainFlank.insert(mainFlank.end(), resAtoms.begin(), resAtoms.end());
    }
    
    for (Residue *candidate: adjacencies.at(residue)) {
        
        Residue *realCandidate;
        if (adjSameResidues) {
            // The candidate is analogous to residue
            realCandidate = candidate;
        } else {
            // The candidate is subsequent to residue
            int candIdx = candidate->getResidueIndexInChain();
            if (candIdx < 1)
                continue;
            realCandidate = candidate->getChain()->getResidues()[candIdx - 1];
        }
        
        // Get the flanking regions of the candidate residue
        vector<Residue *> candChainResidues = realCandidate->getChain()->getResidues();
        int candIdx = realCandidate->getResidueIndexInChain();
        int candChainSize = realCandidate->getChain()->residueSize();
        int candFlankLeft = candIdx >= flankLength ? flankLength : candIdx;
        int candFlankRight = candIdx < candChainSize - flankLength ? flankLength : chainSize - candIdx - 1;
        if (candFlankLeft != realFlankLeft || candFlankRight != realFlankRight)
            continue;
        vector<Atom *> candFlank;
        for (int i = candIdx - candFlankLeft; i < candIdx + candFlankRight; i++) {
            vector<Atom *> resAtoms = candChainResidues[i]->getAtoms();
            candFlank.insert(candFlank.end(), resAtoms.begin(), resAtoms.end());
        }
        
        // Compare RMSD
        float rmsd = RMSDCalculator::rmsd(mainFlank, candFlank);
        if (rmsd <= rmsdCutoff) {
            results.push_back(realCandidate);
        }
    }
    
    return results;
}

void SeedGraph::computeRepresentativeResidues() {
    for (Residue *res: getResidues()) {
        if (_representativeResidues.count(res) != 0)
            continue;
        
        // Set this residue as the representative residue
        for (Residue *relatedRes: neighborhood(res).getResidues()) {
            _representativeResidues[relatedRes] = res;
        }
        _representativeResidues[res] = res;
    }
}

vector<Residue *> SeedGraph::topologicalOrder() {
    deque<Residue *> result;
    unordered_set<Residue *> temporary;
    stack<pair<bool, Residue *>> mStack;
    unordered_set<Residue *> visited;
    
    // Add all residues to the stack preemptively
    for (auto item: adjacencies) {
        Residue *res = item.first;
        mStack.push(make_pair(false, res));
    }

    while (mStack.size() > 0) {
        if (visited.size() % 100 == 0)
            cout << visited.size() << " residues seen so far" << endl;
        
        pair<bool, Residue *> item = mStack.top();
        mStack.pop();
        Residue *res = item.second;
        if (item.first) {
            // All descendants have been processed; add this residue to result
            result.push_front(res);
            temporary.erase(res);
            continue;
        }
        
        if (visited.count(res) != 0)
            continue;
        
        visited.insert(res);
        mStack.push(make_pair(true, res));
        temporary.insert(res);
        
        for (Residue *neighbor: reverseAdjacencies[res]) {
            if (temporary.count(neighbor) != 0) {
                cout << "Found a cycle!?" << endl;
                continue;
            }
            
            mStack.push(make_pair(false, neighbor));
        }
    }
    
    vector<Residue *> retVal;
    retVal.insert(retVal.end(), result.begin(), result.end());
    return retVal;
}
