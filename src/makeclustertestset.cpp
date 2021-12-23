#include "makeclustertestset.h"
#include <algorithm>
#include <random>
#include <math.h>

void clusterTestSet::selectStructures(int negType) {
    vector<string> structureNames = MstUtils::fileToArray(structuresList);
    int truePosCount = 0;
    int trueNegCount = 0;
    int structureCount = 0;
    vector<pair<string, mstreal>> positiveScores;
    vector<pair<string, mstreal>> negativeScores;
    cout << "loading structure cache. List: " << structuresList << " Prefix: " << structuresPath << endl;
    StructureCache* structuresCache = new StructureCache(structuresPath);
    structuresCache->preloadFromPDBList(structuresList);
    auto it = structuresCache->begin();
    ConFind confindTarget(rl, *target);
    int combinedStructureCount = 0;
    while (it != structuresCache->end()) {
        Structure *curStructure = *it;
        string curName = curStructure->getName();
        curName = curName.substr(curName.rfind("/") + 1);
        curName = curName.substr(0, curName.rfind("."));
        mstreal dist = 0;
        int count = 0;
        for (Residue* seedRes : curStructure->getResidues()) {
            Atom* seedCA = seedRes->findAtom("CA");
            mstreal minDist = 10000;
            for (Residue* peptideRes : peptide->getResidues()) {
                Atom* peptideCA = peptideRes->findAtom("CA");
                mstreal curDist = sqrt(pow(seedCA->getX() - peptideCA->getX(), 2) + pow(seedCA->getY() - peptideCA->getY(), 2) + pow(seedCA->getZ() - peptideCA->getZ(), 2));
                if (curDist < minDist)
                    minDist = curDist;  
            }
            dist += minDist;
        }
        dist /= curStructure->getResidues().size();
        pair<string, mstreal> curPair;
        curPair = make_pair(curName, dist);
        positiveScores.emplace_back(curPair);    
        truePosCount++;

        // Get statistics on negatives    
        if ((negType > 0) && (prepareCombinedStructure(curStructure))) {
            combinedStructureCount++;
            if (MstUtils::randInt(structuresCache->size()) <= nStructures) {
                contactParams cParams;
                Chain* poseChain = &((*target)[target->chainSize() - 1]);
                vector<Residue*> poseResidues = poseChain->getResidues();
                
                contactList bbConts; contactList sbConts; contactList bsConts; contactList ssConts;
                // note from craig: switched bsConts and sbConts here as we want sb to be SC from target
                // normal order bbConts, bsConts, sbConts, ssConts
                splitContacts(*target, poseResidues, rl, cParams, false, bbConts, sbConts, bsConts, ssConts);
                contactList conts = contactListUnion({bbConts, sbConts, bsConts, ssConts});
                vector<Residue*> targetRes = conts.srcResidues();

                if (negType == 1) {
                    mstreal avFreedom = 0;
                    for (Residue* res : targetRes) {
                        avFreedom += confindTarget.getFreedom(getTargetResidue(res));
                    }
                    avFreedom /= targetRes.size();
                    pair<string, mstreal> curPair;
                    curPair = make_pair(curName, avFreedom);
                    negativeScores.emplace_back(curPair);   
                } else if (negType == 2) {
                    pair<string, mstreal> curPair;
                    curPair = make_pair(curName, conts.size());
                    negativeScores.emplace_back(curPair);
                }
            }
            resetCombinedStructure();
        }
        structureCount++;
        it++;
    }
    cout << "Number of combined structures: " << combinedStructureCount << "." << endl;

    delete structuresCache;

    // Positive samples
    sort(positiveScores.begin(), positiveScores.end(), [=](std::pair<string, mstreal>& a, std::pair<string, mstreal>& b)
        {
            return a.second < b.second;
        }
    );

    vector<pair<string, mstreal>> possiblePositives;
    for (pair<string, mstreal> structureDistPair : positiveScores) {
        cout << structureDistPair.first << " " << structureDistPair.second << endl;
        if (structureDistPair.second > 1)
            break;
        possiblePositives.emplace_back(structureDistPair);
    }

    distLimit = 0;

    while (truePositives.size() < nStructures/3) {
        bool newStructure = false;
        while (!newStructure) {
            int randInt = MstUtils::randInt(0,possiblePositives.size()-1);
            pair<string, mstreal> structureDistPair = possiblePositives[randInt];
            if (find(truePositives.begin(), truePositives.end(), structureDistPair.first) == truePositives.end()) {
                newStructure = true;
                truePositives.insert(structureDistPair.first);
                if (structureDistPair.second > distLimit)
                    distLimit = structureDistPair.second;
            }
        }
    }

    // Freedom or contact selection
    sort(negativeScores.begin(), negativeScores.end(), [=](std::pair<string, mstreal>& a, std::pair<string, mstreal>& b)
        {
            return a.second < b.second;
        }
    );
    for (pair<string, mstreal> structureFreedomPair : negativeScores) {
        if (find(truePositives.begin(), truePositives.end(), structureFreedomPair.first) == truePositives.end()) 
            trueNegatives.insert(structureFreedomPair.first);
        if (truePositives.size() + trueNegatives.size() >= nStructures*2/3)
            break;
    }

    // Random selection
    int nRemaining = nStructures - truePositives.size() - trueNegatives.size();
    while (unkownSamples.size() < nRemaining) {
        bool newStructure = false;
        while (!newStructure) {
            int randInt = MstUtils::randInt(0,structureNames.size()-1);
            if ((find(truePositives.begin(), truePositives.end(), structureNames[randInt]) == truePositives.end()) && (find(trueNegatives.begin(), trueNegatives.end(), structureNames[randInt]) == trueNegatives.end())) {
                newStructure = true;
                unkownSamples.insert(structureNames[randInt]);
            }
        }
    }
}

bool clusterTestSet::prepareCombinedStructure(Structure *seed) {
    
    AtomPointerVector psTargetAPV = target->getAtoms();
    ProximitySearch ps = ProximitySearch(psTargetAPV, 40.0);
    double vdwRadius = 0.3;

    // Clean up the seed
    Structure tmpStructure(*seed);
    Structure pose;
    cleanStructure(tmpStructure, pose, true, true, true);
    
    if (pose.residueSize() == 0) {
        return false;
    } else if (pose.chainSize() > 1) {
        return false;
    }
    
    removeSideChains(pose); // may not always want to do this
    AtomPointerVector poseAPV = pose.getAtoms();
    set<int> structResis; // empty since we are not excluding any chains
    int adjustNum = 0;
    
    //cout << "Checking for clash..." << endl;
    if (isClash(ps, psTargetAPV, poseAPV, vdwRadius)) {
        // cout << "Clash detected for seed " << seed->getName() << "." << endl;
        return false;
    }
    //cout << "Checked clash" << endl;
    
    // Add the seed to the target as a new chain
    Chain* c = new Chain(generateChainID(target->chainSize()), "", pose.getResidues());
    for (int j = 0; j < c->residueSize(); j++) {
        Residue* res = &c->getResidue(j);
        res->setNum(j + 1);
    }
    target->appendChain(c); // add chain c to the target
    
    return true;
}

void clusterTestSet::resetCombinedStructure() {
    Chain* poseChain = &((*target)[target->chainSize() - 1]);
    target->deleteChain(poseChain);
}

void clusterTestSet::writeStructures(string& outputPath) {
    if (MstSys::fileExists(outputPath)) MstSys::crm(outputPath);
    fstream out;
    MstUtils::openFile(out, outputPath, ios_base::out);
    for (string structureName : truePositives) {
        out << structureName << endl;
    }
    out << endl;
    for (string structureName : unkownSamples) {
        out << structureName << endl;
    }
    out << endl;
    for (string structureName : trueNegatives) {
        out << structureName << endl;
    }
    out.close();
}