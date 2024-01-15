#include "pathsampler.h"
#include "secondarystructure.h"
#include <chrono>
#include <list>

using namespace std::chrono;

void printPath(vector<Residue *> &path) {
    for (int i = 0; i < path.size(); i++) {
        cout << path[i]->getStructure()->getName() << ":" << path[i]->getResidueIndexInChain();
        if (i < path.size() - 1)
            cout << ",";
        else
            cout << endl;
    }
}


// ============= PathSampler (parent class) =============


bool PathSampler::emplacePathFromResidues(vector<Residue *> path, vector<PathResult> &results, set<Residue*> fixedResidues, bool ignore_clashes) {
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

    Structure fused; fusionOutput fuserScore; vector<vector<Residue*>> topologySeedResidues;
    int seedStartIdx = fusePath(path, fused, fuserScore, topologySeedResidues, fixedResidues);
    int interchain_clash = 0; int intrachain_clash = 0;
    bool has_clash = pathClashes(fused, seedStartIdx, interchain_clash, intrachain_clash);
    
    if (has_clash && !ignore_clashes) return false;

    results.emplace_back(path, fused, topologySeedResidues, seedStartIdx, fuserScore, interchain_clash, intrachain_clash);
    
    if (has_clash) return false;
    else return true;
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

vector<Residue *> PathSampler::getTargetResidues(const vector<Residue *> &pathResidues) {
//    unordered_set<pair<string, int>, pair_hash> targetResidueCodes;
    /**
     When writing seed + anchor (extended fragment) files the residue number of the anchor residues were set to the residue index in the structure
    */
    set<int> targetResidueIdx;
    Structure *lastStructure = nullptr;
    for (Residue *res: pathResidues) {
        if (res->getStructure() == lastStructure)
            continue;
        Structure *s = res->getStructure();
        for (int i = 0; i < s->chainSize(); i++) {
            Chain &chain = s->getChain(i);
            if (chain.getID() == "0") continue;
            for (int j = 0; j < chain.residueSize(); j++) {
                targetResidueIdx.insert(chain.getResidue(j).getNum());
            }
        }
        lastStructure = s;
    }

    // Sort and get from target
    vector<Residue *> targetResidues;
    for (int i = 0; i < _target->residueSize(); i++) {
        Residue &res = _target->getResidue(i);
        if (targetResidueIdx.count(res.getResidueIndex()) != 0) {
            targetResidues.push_back(&res);
        }
    }
    return targetResidues;
}

pair<vector<Residue *>, vector<int>> PathSampler::getMappedMatchResidues(const Structure &seedStructure, const map<int,int> &targetPositions) {
    vector<Residue *> residues;
    vector<int> indexes;

    for (int i = 0; i < seedStructure.chainSize(); i++) {
        Chain &chain = seedStructure.getChain(i);
        if (chain.getID() == "0") continue;
        for (int j = 0; j < chain.residueSize(); j++) {
            Residue* R = &chain.getResidue(j);
            residues.push_back(R);
            MstUtils::assertCond(targetPositions.count(R->getNum()) != 0, "Missing target position with index: " + to_string(R->getNum()),"PathSampler::getMappedMatchResidues");
            indexes.push_back(targetPositions.at(R->getNum()));
        }
    }

    return make_pair(residues, indexes);
}

int PathSampler::fusePath(const vector<Residue *> &residues, Structure &fusedPath, fusionOutput& scores, vector<vector<Residue*>>& topologySeedResidues, set<Residue*> fixedResidues) {
    vector<Residue *> targetResidues = getTargetResidues(residues);
    fusionTopology topology(residues.size() + targetResidues.size());

    // Add target residues (there are anchor positions)
    int topologyIndex = 0;
    map<int,int> targetPositions;
    vector<int> targetIndexes;
    for (int i = 0; i < targetResidues.size(); i++) {
        targetIndexes.push_back(i);
        Residue *res = targetResidues[i];
        // map residue index in target structure to index in topology
        targetPositions[res->getResidueIndex()] = i;
    }
    if (!targetResidues.empty()) {
        topology.addFragment(targetResidues, targetIndexes);
        topology.addFixedPositions(targetIndexes);
        cout << "Added " << targetResidues.size() << " target residues, " << targetIndexes.front() << " to " << targetIndexes.back() << endl;
        topologyIndex += targetResidues.size();
    }
    
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
        
        //if residues from this seed are in fixed residues, up-weight the whole fragment
        mstreal weight = 1.0;
        if (fixedResidues.count(allResidues[allResidues.size()-1]) > 0) {
            weight = 100000.0;
            cout << "add residues from the peptide with higher weight" << endl;
        }
        
        topology.addFragment(allResidues, allIndexes, weight);
        topologySeedResidues.push_back(allResidues);
        topologyIndex += seedRes.size();
    }
    if (twoStepFuse) {
        cout << "two step fuse" << endl;
        fusionParams opts;
        vector<mstreal> ics = opts.getIntCoorFCs();
        opts.setIntCoorFCs({0, 0, 0});
        Structure initialFusedPath = Fuser::fuse(topology, scores, opts);
        opts.setIntCoorFCs(ics);
        opts.setStartingStructure(initialFusedPath);
        fusedPath = Fuser::fuse(topology, scores, opts);
    } else fusedPath = Fuser::fuse(topology, scores);
    return seedStartIdx;
}

bool PathSampler::pathClashes(const Structure &path, int seedStartIdx, int &interchain_clash, int &intrachain_clash) {
    Structure pathOnly;
    Chain *c = new Chain;
    for (int i = seedStartIdx; i < path.residueSize(); i++) {
        c->appendResidue(new Residue(path.getResidue(i)));
    }
    pathOnly.appendChain(c);

    AtomPointerVector pathAPV(pathOnly.getAtoms());
    
    set<pair<Residue*, Residue*> > exclude;
    interchain_clash = numClash(ps, targetAPV, pathAPV, exclude, 0.7, -1);
    if (interchain_clash > 0) cout << "Path clashes with protein" << endl;
    
    ProximitySearch pathPS(pathAPV, 10.0);
    exclude.clear();
    intrachain_clash = numClash(pathPS, pathAPV, pathAPV, exclude, 0.7, -1);
    if (intrachain_clash > 0) cout << "Path clashes with self" << endl;
    if (interchain_clash > 0 || intrachain_clash > 0) return true;
    else return false;
}


// ============= SeedGraphPathSampler =============


vector<PathResult> SeedGraphPathSampler::sample(int numPaths) {
    vector<PathResult> results;

    while (results.size() < numPaths) {
        Residue *startRes;

        // the position of the "first" residue of the seed that was most recently added to the path
        // this residue represents the boundary: if this residue cannot be extended (and thus all
        // before it could not be either), then the path cannot continue growing in the current direction
        int lastSeedLimit;
        vector<Residue*> initial_path;
        vector<Residue*> path;
        if (!_initialPathResidues.empty()) {
            // The path begins with the already selected residues
            startRes = _initialPathResidues.back();
            path = startRes->getParent()->getResidues();
            // the path *must* include these residues, so the C-terminal limit is set by this
            auto it = find(path.begin(),path.end(),startRes);
            lastSeedLimit = it - path.begin();
        } else {
            // Pick a residue at random from the graph
            startRes = _startingResidues[MstUtils::randInt(_startingResidues.size())];
            path = startRes->getParent()->getResidues();
            lastSeedLimit = 0;
        }
        string seedName = startRes->getStructure()->getName();
        int originalCount = path.size();
        int firstSeedEnd = -1; //the position of the last residue of the first seed that was added to the path
        initial_path = path;
        
        // Extend the C-terminus of the path
        int endIndex = path.size() - 1;
        vector<Residue *> adj;
        while (endIndex >= lastSeedLimit) {
            unordered_set<Residue *> allAdj = _graph->forwardNeighbors(path[endIndex]);

            // ensure that the adjacent residues are part of new seeds
            adj.clear();
            for (Residue *res: allAdj) {
                if (res->getStructure()->getName() == seedName) continue;
                adj.push_back(res);
            }

            if (!adj.empty()) {
                // the are potential residue(s) on other seeds that can be added to the path, select
                // one seed to add
                Residue *nextRes = adj[MstUtils::randInt(adj.size())];
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
                // there were no residues on a different seed to add to the path, check if there are
                // potential bonds with the neighbouring residue in the N-terminal direction
                endIndex--;
            }
        }

        // Extend the N-terminus of the path
        endIndex = 0;
        if (!_initialPathResidues.empty()) {
            // the last seed limit is the position of the first required residue in the path
            Residue* first = _initialPathResidues.front();
            auto it = find(path.begin(),path.end(),first);
            lastSeedLimit = it - path.begin();
        } else {
            // the last seed limit is either the position of the last residue from the first seed
            // or just the last residue of the first seed, if it was never extended towards the C-terminus
            lastSeedLimit = (firstSeedEnd >= 0 ? firstSeedEnd : originalCount - 1);
        }
        seedName = path[0]->getStructure()->getName();
        adj.clear();
        while (endIndex <= lastSeedLimit) {
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
                lastSeedLimit = resIdx - newResidues.begin();
                endIndex = 0;
            } else {
                endIndex++;
            }
        }
        
        attempts++;
        
        // Check if path is accepted
        if (!acceptSingleSeedPaths && (path == initial_path)) {
            cout << "Sampled path is the same as the initially selected seed" << endl;
            no_overlaps++;
            continue;
        } else if (_sampledPaths.find(path) != _sampledPaths.end()){
            cout << "Path is redundant to a previously sampled one " << endl;
            redundant++;
            continue;
        } else if ((path.size() < minimumLength)) {
            cout << "Sampled path is too short: " << path.size() << endl;
            too_short++;
            continue;
        } else if (!emplacePathFromResidues(path, results)) {
            cout << "Sampled path clashes" << endl;
            clashes++;
            continue;
        }
        
        _sampledPaths.insert(path);

        if (results.size() % 100 == 0) {
            cout << results.size() << " paths generated" << endl;
        }
    }
    return results;
}

vector<PathResult> SeedGraphPathSampler::fusePaths(const vector<string>& path_specifiers) {
    vector<PathResult> results;
    
    int path_count = 0;
    for (string path_spec : path_specifiers) {
        path_count++;
        vector<Residue*> path = pathResiduesFromSpecifier(path_spec);
        if (path.empty()) MstUtils::error("Path has no residues");
        bool ignore_clashes = true;
        bool no_clash = emplacePathFromResidues(path, results, _fixedResidues, ignore_clashes);
        if (!no_clash) cout << "Path " << path_count << " has a clash (or is missing atoms)" << endl;
    }
    return results;
}

vector<Residue*> SeedGraphPathSampler::pathResiduesFromSpecifier(string path_spec) {
    vector<Residue*> path;
    //split name into SeedGraph residue_names
    vector<string> residue_names = MstUtils::split(path_spec,";");
    for (string residue_name : residue_names) {
        Residue* R = _graph->getResidueFromFile(residue_name);
        if (R == nullptr) MstUtils::error("Residue not found in seed binary file: "+residue_name,"SeedGraphPathSampler::pathResiduesFromSpecifier");
        path.push_back(R);
    }
    return path;
}


// ============= ClusterTreePathSampler =============


vector<PathResult> ClusterTreePathSampler::sample(int numPaths) {
    vector<PathResult> results;
    
    // Prevents the search tree from searching for unnecessary results.
    // Most seed residues will have much fewer than 100 matches.
    _searchTree->sufficientMatches = 100;
    while (results.size() < numPaths) {
        // Pick a residue at random
        Residue *startRes = _searchTreeResidues[MstUtils::randInt(_searchTreeResidues.size())];
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

        if ((path.size() <= originalCount) || (path.size() < minimumLength)) {
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

    return results;
}

mstreal ClusterTreePathSampler::cosineAngle(const vector<Residue *> &res1, const vector<Residue *> &res2) {
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

vector<int> ClusterTreePathSampler::shuffleResultIndexes(ClusterSearchResults &searchResults) {
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

Residue * ClusterTreePathSampler::findSeedToAdd(const vector<Residue *> &path, int endIndex, const string &currentSeedName, bool appendingForward) {
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

