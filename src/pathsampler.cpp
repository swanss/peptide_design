#include "pathsampler.h"
#include "Util.h"
#include "secondarystructure.h"
#include <chrono>
#include <list>

using namespace std::chrono;

vector<PathResult> PathSampler::sample(int numPaths) {
    vector<PathResult> paths;
    if (_graph != nullptr)
        sampleFromGraph(numPaths, paths);
    else if (_searchTree != nullptr)
        sampleFromOverlapTree(numPaths, paths);
    return paths;
}

void PathSampler::sampleFromGraph(int numPaths, vector<PathResult> &results) {
    // This does not support new features like unique seed restrictions
    unordered_set<Residue *> residues = _graph->getResidues();
    vector<Residue *> residuesOrdered(residues.begin(), residues.end());

    while (results.size() < numPaths) {
        // Pick a residue at random from the graph
        Residue *startRes = residuesOrdered[rand() % residuesOrdered.size()];
        vector<Residue *> path(startRes->getParent()->getResidues());
        string seedName = startRes->getStructure()->getName();
        int originalCount = path.size();
        int firstSeedEnd = -1;

        // Extend forward
        int endIndex = path.size() - 1;
        int lastSeedLimit = 0;
        vector<Residue *> adj;
        while (endIndex > lastSeedLimit || !adj.empty()) {
            unordered_set<Residue *> allAdj = _graph->forwardNeighbors(path[endIndex]);
            adj.clear();
            for (Residue *res: allAdj) {
                if (res->getStructure()->getName() == seedName) continue;
                adj.push_back(res);
            }

            if (!adj.empty()) {
                Residue *nextRes = adj[rand() % adj.size()];
                seedName = nextRes->getStructure()->getName();
                path.erase(path.begin() + endIndex + 1, path.end());
                vector<Residue *> newResidues = nextRes->getParent()->getResidues();
                auto resIdx = find(newResidues.begin(), newResidues.end(), nextRes);
                if (resIdx == newResidues.end())
                    MstUtils::error("Residue not found in its own chain!");
                path.insert(path.end(), resIdx, newResidues.end());
                if (firstSeedEnd < 0)
                    firstSeedEnd = endIndex;
                lastSeedLimit = endIndex + 1;
                endIndex = path.size() - 1;
            } else {
                endIndex--;
            }
        }
    
        // Extend backward
        endIndex = 0;
        lastSeedLimit = (firstSeedEnd >= 0 ? firstSeedEnd : originalCount); 
        seedName = path[0]->getStructure()->getName();
        adj.clear();
        while (endIndex < lastSeedLimit || !adj.empty()) {
            unordered_set<Residue *> allAdj = _graph->backwardNeighbors(path[endIndex]);
            adj.clear();
            for (Residue *res: allAdj) {
                if (res->getStructure()->getName() == seedName) continue;
                adj.push_back(res);
            }

            if (!adj.empty()) {
                Residue *nextRes = adj[rand() % adj.size()];
                seedName = nextRes->getStructure()->getName();
                path.erase(path.begin(), path.begin() + endIndex);
                vector<Residue *> newResidues = nextRes->getParent()->getResidues();
                auto resIdx = find(newResidues.begin(), newResidues.end(), nextRes);
                if (resIdx == newResidues.end())
                    MstUtils::error("Residue not found in its own chain!");
                path.insert(path.begin(), newResidues.begin(), resIdx + 1);
                lastSeedLimit = distance(newResidues.begin(), resIdx) + 1;
                endIndex = 0;
            } else {
                endIndex++;
            }
        }

        if (path.size() <= originalCount)
            continue;
        if (!emplacePathFromResidues(path, results))
            continue;
        if (results.size() % 100 == 0)
            cout << results.size() << " paths generated" << endl;
    }
}

AtomPointerVector PathSampler::buildAPV(const vector<Residue *>::const_iterator &begin, const vector<Residue *>::const_iterator &end) {
    AtomPointerVector result;
    for (auto it = begin; it != end; ++it) {
        Residue *res = *it;
        for (Atom *a: res->getAtoms()) {
            result.push_back(a);
        }
    }
    return result;
}

void printPath(vector<Residue *> &path) {
    for (int i = 0; i < path.size(); i++) {
        cout << path[i]->getStructure()->getName() << ":" << path[i]->getResidueIndexInChain();
        if (i < path.size() - 1)
            cout << ",";
        else
            cout << endl;
    }
}

mstreal PathSampler::cosineAngle(const vector<Residue *> &res1, const vector<Residue *> &res2) {
    Atom *first1 = res1.front()->findAtom("CA");
    Atom *last1 = res1.back()->findAtom("CA");
    Atom *first2 = res2.front()->findAtom("CA");
    Atom *last2 = res2.back()->findAtom("CA");
    mstreal dx1 = last1->getX() - first1->getX();
    mstreal dx2 = last2->getX() - first2->getX();
    mstreal dy1 = last1->getY() - first1->getY();
    mstreal dy2 = last2->getY() - first2->getY();
    mstreal dz1 = last1->getZ() - first1->getZ();
    mstreal dz2 = last2->getZ() - first2->getZ();
    mstreal dot = dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
    mstreal mag1 = sqrt(pow(dx1, 2.0) + pow(dy1, 2.0) + pow(dz1, 2.0));
    mstreal mag2 = sqrt(pow(dx2, 2.0) + pow(dy2, 2.0) + pow(dz2, 2.0));
    return dot / (mag1 * mag2);
}

Residue *PathSampler::findSeedToAdd(const vector<Residue *> &path, int endIndex, const string &currentSeedName, bool appendingForward) {
    // Build a template on which we want to find overlaps
    vector<Residue *> overlapSource;
    if (appendingForward)
        overlapSource.insert(overlapSource.end(), path.begin() + endIndex + 1 - overlapLength, path.begin() + endIndex + 1);
    else
        overlapSource.insert(overlapSource.end(), path.begin() + endIndex, path.begin() + endIndex + overlapLength);

    AtomPointerVector apv = buildAPV(overlapSource.begin(), overlapSource.end());
    ClusterSearchResults searchResults = _searchTree->search(apv, overlapRMSD);

    if (searchResults.size() == 0)
        return nullptr;

    Residue *currentRes = nullptr;
    vector<Residue *> overlapResidues;

    // Shuffle indexes to try
    vector<int> indexes = shuffleResultIndexes(searchResults);

    int resultIndex = 0;
    while (resultIndex < searchResults.size()) {
        overlapResidues = searchResults.getResultResidues(indexes[resultIndex++]);
        Residue *res = appendingForward ? overlapResidues.back() : overlapResidues.front();
        string newSeedName = res->getStructure()->getName();

        // Check that this residue is not part of the current seed
        if (newSeedName == currentSeedName)
            continue;
        
        // If uniqueSeeds is enabled, check that it is a unique seed
        if (uniqueSeeds && _usedSeeds.count(newSeedName) != 0)
            continue;

        // Make sure this residue has enough residues to add
        Chain *seedChain = res->getParent();
        if (appendingForward && endIndex + seedChain->residueSize() - res->getResidueIndexInChain() <= (int)path.size())
            continue;
        else if (!appendingForward && res->getResidueIndexInChain() - endIndex <= 0)
            continue;

        // Make sure the dot product of the overlap is sufficiently high
        if (cosineAngle(overlapSource, overlapResidues) < minCosineAngle)
            continue;
        currentRes = res;
        break;
    }

    return currentRes;
}

vector<int> PathSampler::shuffleResultIndexes(ClusterSearchResults &searchResults) {
    vector<int> indexes(searchResults.size());

    if (preferredSecondaryStructure != nullptr) {
        // Weighted random
        secondaryStructureClassifier classifier(7, 8);
        list<pair<int, float>> weights;
        float totalWeight = 0.0f;
        for (int i = 0; i < searchResults.size(); i++) {
            Structure *s = searchResults.getFullStructure(i);
            // Upweight if the correct ss type is a majority of the seed
            int correctSSCount = 0;
            for (Residue *res: s->getResidues()) {
                if (classifier.classifyResidue(res) == *preferredSecondaryStructure) {
                    correctSSCount++;
                }    
            }
            float weight = correctSSCount >= (s->residueSize() / 2) ? secondaryStructureWeight : 1.0f;
            totalWeight += weight;
            weights.push_back(make_pair(i, weight));
        }

        for (int i = 0; i < indexes.size(); i++) {
            float val = rand() / (RAND_MAX / totalWeight);
            float currThreshold = 0.0f;
            for (auto it = weights.begin(); it != weights.end(); ++it) {
                float weight = (*it).second;
                currThreshold += weight;
                if (val <= currThreshold) {
                    indexes[i] = (*it).first;
                    weights.erase(it);
                    totalWeight -= weight;
                    break;
                }
            }
        }

    } else {
        // Uniform random
        iota(indexes.begin(), indexes.end(), 0);
        random_shuffle(indexes.begin(), indexes.end());
    }
    return indexes;
}

void PathSampler::sampleFromOverlapTree(int numPaths, vector<PathResult> &results) {
    // Prevents the search tree from searching for unnecessary results.
    // Most seed residues will have much fewer than 100 matches.
    _searchTree->sufficientMatches = 100;
    while (results.size() < numPaths) {
        // Pick a residue at random
        Residue *startRes = _searchTreeResidues[rand() % _searchTreeResidues.size()];
        vector<Residue *> path(startRes->getParent()->getResidues());
        string seedName = startRes->getStructure()->getName();
        if (uniqueSeeds && _usedSeeds.count(seedName) != 0) {
            cout << "Seed already used" << endl;
            continue;
        }

        int originalCount = path.size();
        //cout << "Original path: ";
        //printPath(path);
        int firstSeedEnd = -1;

        // Extend forward
        int endIndex = path.size() - 1;
        int lastSeedLimit = max(0, overlapLength - 2);
        bool hasResults = false;
        while (endIndex > lastSeedLimit || hasResults) {
            Residue *currentRes = findSeedToAdd(path, endIndex, seedName, true);
            if (currentRes == nullptr) {
                // None of the results work
                hasResults = false;   
                endIndex--;
                continue;
            }

            // Add the end of currentRes's seed
            path.erase(path.begin() + endIndex + 1, path.end());
            seedName = currentRes->getStructure()->getName();
            vector<Residue *> seedResidues = currentRes->getParent()->getResidues();
            path.insert(path.end(), seedResidues.begin() + currentRes->getResidueIndexInChain() + 1, seedResidues.end());

            if (firstSeedEnd < 0)
                firstSeedEnd = endIndex;
            //cout << "Forward (" << endIndex << "): ";
            //printPath(path);
            lastSeedLimit = endIndex + 1;
            endIndex = path.size() - 1;
        }
    
        // Extend backward
        endIndex = 0;
        lastSeedLimit = (firstSeedEnd >= 0 ? firstSeedEnd : originalCount) - overlapLength + 1; 
        seedName = path[0]->getStructure()->getName();
        hasResults = false;
        while (endIndex < lastSeedLimit || hasResults) {
            Residue *currentRes = findSeedToAdd(path, endIndex, seedName, false);
            if (currentRes == nullptr) {
                // None of the results work
                hasResults = false;   
                endIndex++;
                continue;
            }
            
            // Add the beginning of currentRes's seed
            path.erase(path.begin(), path.begin() + endIndex);
            seedName = currentRes->getStructure()->getName();
            vector<Residue *> seedResidues = currentRes->getParent()->getResidues();
            path.insert(path.begin(), seedResidues.begin(), seedResidues.begin() + currentRes->getResidueIndexInChain());

            //cout << "Backward (" << endIndex << "): ";
            //printPath(path);
            lastSeedLimit = currentRes->getResidueIndexInChain();
            endIndex = 0;
        }

        if (path.size() <= originalCount) {
            cout << "Path too small" << endl;
            continue;
        }
        if (!emplacePathFromResidues(path, results)) {
            cout << "Path clashes" << endl;
            continue;
        }
        //if (results.size() % 100 == 0)
        if (uniqueSeeds) {
            for (Residue *res: path)
                _usedSeeds.insert(res->getStructure()->getName());

            cout << results.size() << " paths generated, " << _usedSeeds.size() << " seeds used" << endl;
        } else
            cout << results.size() << " paths generated" << endl;
    }
}

bool PathSampler::emplacePathFromResidues(vector<Residue *> path, vector<PathResult> &results) {
    // From mstfuser.cpp: we require that every residue in the path has these four atoms
    static vector<string> bba = {"N", "CA", "C", "O"};

    bool isValid = true;
    for (Residue *res: path) {
        for (string atomCode: bba) {
            if (!res->atomExists(atomCode)) {
                isValid = false;
                break;
            }
        }
        if (!isValid)
            break;
    }
    if (!isValid)
        return false;

    Structure fused;
    int seedStartIdx = fusePath(path, fused);
    if (pathClashes(fused, seedStartIdx))
        return false;

    results.emplace_back(path, fused, seedStartIdx);
    return true;
}

bool PathSampler::pathClashes(const Structure &path, int seedStartIdx) {
    Structure pathOnly;
    Chain *c = new Chain;
    for (int i = seedStartIdx; i < path.residueSize(); i++) {
        c->appendResidue(new Residue(path.getResidue(i)));
    }
    pathOnly.appendChain(c);

    AtomPointerVector pathAPV(pathOnly.getAtoms());
    if (isClash(ps, targetAPV, pathAPV, 0.7)) {
        return true;
    }
    ProximitySearch pathPS(pathAPV, 10.0);
  //not clear where this 
//    if (isClashSingleStructure(pathPS, pathAPV, 0.7)) {
//        return true;
//    }
    return false;
}

vector<Residue *> PathSampler::getTargetResidues(const vector<Residue *> &residues) {
    unordered_set<pair<string, int>, pair_hash> targetResidueCodes;
    Structure *lastStructure = nullptr;
    for (Residue *res: residues) {
        if (res->getStructure() == lastStructure)
            continue;
        Structure *s = res->getStructure();
        for (int i = 0; i < s->chainSize(); i++) {
            Chain &chain = s->getChain(i);
            if (chain.getID() == "0") continue;
            for (int j = 0; j < chain.residueSize(); j++) {
                // Refer to target residue by its labeled number in the match, NOT residue index
                targetResidueCodes.insert(make_pair(chain.getID(), chain.getResidue(j).getNum()));
            }
        }
        lastStructure = s;
    } 

    // Sort and get from target
    vector<Residue *> targetResidues;
    for (int i = 0; i < _target->residueSize(); i++) {
        Residue &res = _target->getResidue(i);
        // Retrieve residue from target by index in chain, NOT number
        if (targetResidueCodes.count(make_pair(res.getChain()->getID(), res.getResidueIndexInChain())) != 0) {
            targetResidues.push_back(&res);
        }
    }
    return targetResidues;
}

pair<vector<Residue *>, vector<int>> PathSampler::getMappedMatchResidues(const Structure &seedStructure, const unordered_map<pair<string, int>, int, pair_hash> &targetPositions) {
    vector<Residue *> residues;
    vector<int> indexes;

    for (int i = 0; i < seedStructure.chainSize(); i++) {
        Chain &chain = seedStructure.getChain(i);
        if (chain.getID() == "0") continue;
        for (int j = 0; j < chain.residueSize(); j++) {
            residues.push_back(&chain.getResidue(j));
            pair<string, int> residueCode(chain.getID(), chain.getResidue(j).getNum());
            MstUtils::assert(targetPositions.count(residueCode) != 0, "Missing target position: " + residueCode.first + to_string(residueCode.second));
            indexes.push_back(targetPositions.at(residueCode));
        }
    }

    return make_pair(residues, indexes);
}

int PathSampler::fusePath(const vector<Residue *> &residues, Structure &fusedPath) {
    vector<Residue *> targetResidues = getTargetResidues(residues);
    fusionTopology topology(residues.size() + targetResidues.size());

    // Add target residues
    int topologyIndex = 0;
    unordered_map<pair<string, int>, int, pair_hash> targetPositions;
    vector<int> targetIndexes;
    for (int i = 0; i < targetResidues.size(); i++) {
        targetIndexes.push_back(i);
        Residue *res = targetResidues[i];
        targetPositions[make_pair(res->getChain()->getID(), res->getResidueIndexInChain())] = i;
    }
    topology.addFragment(targetResidues, targetIndexes);
    topology.addFixedPositions(targetIndexes);
    //cout << "Added " << targetResidues.size() << " target residues, " << targetIndexes.front() << " to " << targetIndexes.back() << endl;
    topologyIndex += targetResidues.size();
    
    // Split path residues by which seed they came from
    vector<vector<Residue *>> residuesBySeed;
    for (Residue *res: residues) {
        if (residuesBySeed.empty() || residuesBySeed.back().back()->getStructure()->getName() != res->getStructure()->getName()) {
            residuesBySeed.push_back({ res });
        } else {
            residuesBySeed.back().push_back(res);
        }
    }

    // Add seed and match residues
    int seedStartIdx = topologyIndex;
    for (vector<Residue *> &seedRes: residuesBySeed) {
        Structure *seedStruct = seedRes.back()->getStructure();
        auto matchRes = getMappedMatchResidues(*seedStruct, targetPositions);
        vector<Residue *> matchResidues = matchRes.first;
        vector<int> matchIndexes = matchRes.second;

        vector<Residue *> seedResidues(seedRes);
        vector<int> seedIndexes(seedResidues.size());
        iota(seedIndexes.begin(), seedIndexes.end(), topologyIndex);

        // Elongate in both directions until either the seed boundary is
        // reached, we run out of seed residues, or we've added overlapLength
        // residues
        int overlapIndex = topologyIndex - 1;
        int overlapNum = 0;
        Residue *currentResidue = seedResidues.front()->previousResidue();
        while (overlapIndex >= seedStartIdx && currentResidue != nullptr && overlapNum < overlapLength) {
            seedResidues.insert(seedResidues.begin(), currentResidue);
            seedIndexes.insert(seedIndexes.begin(), overlapIndex);
            overlapIndex--;
            overlapNum++;
            currentResidue = currentResidue->previousResidue();
        }
        overlapIndex = topologyIndex + seedRes.size();
        overlapNum = 0;
        currentResidue = seedResidues.back()->nextResidue();
        while (overlapIndex < seedStartIdx + residues.size() && currentResidue != nullptr && overlapNum < overlapLength) {
            seedResidues.push_back(currentResidue);
            seedIndexes.push_back(overlapIndex);
            overlapIndex++;
            overlapNum++;
            currentResidue = currentResidue->nextResidue();
        } 

        //cout << "Adding " << seedResidues.size() << " residues from " << seedStruct->getName() << ": " << seedIndexes.front() << " to " << seedIndexes.back() << endl;
        vector<Residue *> allResidues(matchResidues);
        allResidues.insert(allResidues.end(), seedResidues.begin(), seedResidues.end());
        vector<int> allIndexes(matchIndexes);
        allIndexes.insert(allIndexes.end(), seedIndexes.begin(), seedIndexes.end());
        topology.addFragment(allResidues, allIndexes);
        topologyIndex += seedRes.size();
    }
    
    fusionOutput scores;
    fusedPath = Fuser::fuse(topology, scores);
    return seedStartIdx;
}
