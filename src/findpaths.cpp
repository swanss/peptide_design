#include <stdio.h>
#include "findpaths.h"

vector<pair<vector<Residue *>, mstreal>> PathFinder::findPaths(int numPaths) {
    vector<Residue *> orderedResidues = graph->topologicalOrder();
    cout << "Computed topological order" << endl;
    
    //scores->computeRepresentativeResidues();
    //cout << "Computed representative residues" << endl;
    
    // Calculate k best paths for all residues, and take the overall best k
    RankList<Residue *> bestStarts(numPaths, false);
    int resIndex = 0;
    for (Residue *res: orderedResidues) {
        if ((++resIndex) % progressPrintInterval == 0) {
            cout << ceil((resIndex / (float)graph->residueSize()) * 100.0f) << "% of residues evaluated" << endl;
        }

        //Residue *rep = scores->representativeResidue(res);
        if (pathPointers.count(res) != 0)
            continue;
        
        computeBestPaths(res, numPaths);
        for (auto item: pathPointers[res].values()) {
            bestStarts.insert(res, item.second);
        }
    }
    cout << "Found all best starts: " << bestStarts.size() << endl;
    
    // Construct the paths corresponding to the best starts
    vector<pair<vector<Residue *>, mstreal>> result;
    auto values = bestStarts.values();
    int i = 0;
    for (pair<Residue *, mstreal> item: values) {
        result.push_back(make_pair(constructPath(bestStarts, i), item.second));
        i++;
    }
    return result;
}

void PathFinder::computeBestPaths(Residue *residue, int numPaths) {
    /*  In pathPointers, store a rank list specifying the next residues leading
        to the best paths starting from this residue. nullptr corresponds to
        the path ending at this residue.
     */
    if (pathPointers.count(residue) > 0) {
        return;
    }
    if (temporaryMarks.count(residue) > 0) {
        cout << "Found cycle around " << graph->writeCodeForResidue(residue) << endl;
        return;
    }
    
    temporaryMarks.insert(residue);
    RankList<Residue *> bestPaths(numPaths, false);
    mstreal *scorePtr = scores->value(residue);
    if (scorePtr == nullptr) {
        pathPointers[residue] = bestPaths;
        return;
    }
    mstreal score = *scorePtr;
    if (score >= 1e9) {
        pathPointers[residue] = bestPaths;
        return;
    }
    
    bestPaths.insert(nullptr, score);
    int numNeighbors = graph->forwardNeighbors(residue).size();
    if (numNeighbors > 1)
        cout << "Has " << numNeighbors << " neighbors" << endl;
    for (Residue *adjacent: graph->forwardNeighbors(residue)) {
        //Residue *adjacentRep = scores->representativeResidue(adjacent);
        if (pathPointers.count(adjacent) == 0) {
            computeBestPaths(adjacent, numPaths);
        }
        
        for (auto item: pathPointers[adjacent].values()) {
            bestPaths.insert(adjacent, item.second + score);
        }
    }
    
//    cout << "Result for " << graph->writeCodeForResidue(residue) << endl;
//    for (auto item: bestPaths.values()) {
//        cout << "\t" << (item.first != nullptr ? graph->writeCodeForResidue(item.first) : "null") << ": " << item.second << endl;
//    }
    pathPointers[residue] = bestPaths;
    temporaryMarks.erase(residue);
}

vector<Residue *> PathFinder::constructPath(RankList<Residue *> &rankList, int rankIndex) {
    vector<Residue *> path;
    if (rankIndex >= rankList.size()) {
        cout << "Rank index out of bounds" << endl;
        return path;
    }
    Residue *startResidue = rankList[rankIndex];
    if (startResidue == nullptr)
        return path;

    path.push_back(startResidue);
    
    // Find the number of times k that this residue has appeared before now
    int i = 0;
    int k = 0;
    for (Residue *res: rankList.items()) {
        if (i == rankIndex) break;
        if (res == startResidue) {
            k++;
        }
        i++;
    }

    // Get the kth best path through that residue
    vector<Residue *> pathEnd = constructPath(pathPointers[startResidue], k);
    path.insert(path.end(), pathEnd.begin(), pathEnd.end());
    return path;
}
