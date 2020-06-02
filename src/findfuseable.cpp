#include "findfuseable.h"
#include "structure_iter.h"
#include <regex>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <float.h>
#include <chrono>
#include <unordered_set>

using namespace std;
using namespace MST;
using namespace std::chrono;

vector<Structure *> FuseCandidateFinder::findNearbySeeds(Structure *seed, vector<Structure *> &structures, float inset)
{
    vector<Structure *> result;
    
    // Compute the bounding box of this seed
    mstreal xlo, xhi, ylo, yhi, zlo, zhi;
    //ProximitySearch::calculateExtent(*seed, xlo, ylo, zlo, xhi, yhi, zhi);
    Atom &minAtom = seed->getResidue(0).getAtom(0);
    Atom &midAtom = seed->getResidue(seed->residueSize() / 2).getAtom(0);
    Atom &maxAtom = seed->getResidue(seed->residueSize() - 1).getAtom(0);
    xlo = min(minAtom.getX(), min(midAtom.getX(), maxAtom.getX()));
    xhi = max(minAtom.getX(), max(midAtom.getX(), maxAtom.getX()));
    ylo = min(minAtom.getY(), min(midAtom.getY(), maxAtom.getY()));
    yhi = max(minAtom.getY(), max(midAtom.getY(), maxAtom.getY()));
    zlo = min(minAtom.getZ(), min(midAtom.getZ(), maxAtom.getZ()));
    zhi = max(minAtom.getZ(), max(midAtom.getZ(), maxAtom.getZ()));

    // Add the inset
    xlo += inset; ylo += inset; zlo += inset;
    xhi -= inset; yhi -= inset; zhi -= inset;
    
    for (Structure *candidate : structures) {
        if (candidate->getName() == seed->getName())
            continue;
        mstreal cand_xlo, cand_xhi, cand_ylo, cand_yhi, cand_zlo, cand_zhi;
        //ProximitySearch::calculateExtent(*candidate, cand_xlo, cand_ylo, cand_zlo, cand_xhi, cand_yhi, cand_zhi);
        Atom &candMin = candidate->getResidue(0).getAtom(0);
        Atom &candMid = candidate->getResidue(candidate->residueSize() / 2).getAtom(0);
        Atom &candMax = candidate->getResidue(candidate->residueSize() - 1).getAtom(0);
        cand_xlo = min(candMin.getX(), min(candMid.getX(), candMax.getX()));
        cand_xhi = max(candMin.getX(), max(candMid.getX(), candMax.getX()));
        cand_ylo = min(candMin.getY(), min(candMid.getY(), candMax.getY()));
        cand_yhi = max(candMin.getY(), max(candMid.getY(), candMax.getY()));
        cand_zlo = min(candMin.getZ(), min(candMid.getZ(), candMax.getZ()));
        cand_zhi = max(candMin.getZ(), max(candMid.getZ(), candMax.getZ()));

        // If bounding boxes overlap, return this seed
        if (intervalsOverlap(xlo, xhi, cand_xlo, cand_xhi) &&
            intervalsOverlap(ylo, yhi, cand_ylo, cand_yhi) &&
            intervalsOverlap(zlo, zhi, cand_zlo, cand_zhi)) {
            result.push_back(candidate);
        }
    }
    
    return result;
}

vector<mstreal> FuseCandidateFinder::computeBatchBoundingBox(vector<Structure *> seeds, double inset) {
    // Compute the bounding box of this seed
    mstreal xloMin, xhiMax, yloMin, yhiMax, zloMin, zhiMax;

    for (Structure *seed: seeds) {
        mstreal xlo, xhi, ylo, yhi, zlo, zhi;
        ProximitySearch::calculateExtent(*seed, xlo, ylo, zlo, xhi, yhi, zhi);
        // Add the inset
        xlo += inset; ylo += inset; zlo += inset;
        xhi -= inset; yhi -= inset; zhi -= inset;

        xloMin = min(xlo, xloMin);
        xhiMax = max(xhi, xhiMax);
        yloMin = min(ylo, yloMin);
        yhiMax = max(yhi, yhiMax);
        zloMin = min(zlo, zloMin);
        zhiMax = max(zhi, zhiMax);
    }
    
    vector<mstreal> result;
    result.push_back(xloMin);
    result.push_back(xhiMax);
    result.push_back(yloMin);
    result.push_back(yhiMax);
    result.push_back(zloMin);
    result.push_back(zhiMax);
    return result;
}

mstreal FuseCandidateFinder::quickRMSD(const vector<Atom *> &atoms1, const vector<Atom *> &atoms2, mstreal threshold) {
    size_t M = atoms1.size();
    size_t N = atoms2.size();
    if (M != N)
        MstUtils::error("Atom lists are different sizes", "FuseCandidateFinder::rmsdWithinThreshold");

    if (atoms1[0]->distance2(atoms2[0]) > threshold * threshold * N)
        return 1e10;

    mstreal residual = 0.0;
    for (int i = 0; i < N; i++) {
        residual += atoms1[i]->distance2(atoms2[i]);
    }
    return sqrt(residual / N); 
}

vector<FuseCandidate> FuseCandidateFinder::findOverlappingSeeds(Structure *seed, vector<Structure *> &structures, bool filterByBBox)
{
    // First filter by bounding box overlap
    vector<Structure *> fuseCandidates = filterByBBox ? findNearbySeeds(seed, structures, -2.0) : structures;
    
    RMSDCalculator rms;
    vector<FuseCandidate> result;
    
    //cout << fuseCandidates.size() << " candidates for alignment check" << endl;
    for (Structure *fuseCandidate : fuseCandidates) {
        
        for (int i = 0; i < fuseCandidate->chainSize(); i++) {
            Chain candChain = fuseCandidate->getChain(i);
            if (candChain.residueSize() < numResOverlap) continue;
            
            for (int j = 0; j < seed->chainSize(); j++) {
                Chain baseChain = seed->getChain(j);
                if (baseChain.residueSize() < numResOverlap) continue;
                
                // Calculate RMSD of chain i of candidate with chain j of base
                if (searchMode == general) {
                    // Scan overlaps from both sides
                    vector<tuple<int, int, float>> bestOverlaps = scanOverlaps(baseChain, candChain);
                    for (auto bestOverlap: bestOverlaps) { // if (get<2>(bestOverlap) <= rmsdCutoff) {
                        int baseOverlap = get<0>(bestOverlap);
                        int candOverlap = get<1>(bestOverlap);
                        
                        /*
                         There are two possible fusions:
                         base - A \_._/ B
                         cand - C /   \ D
                         A.D (base -> cand) and C.B (cand -> base).
                         If the overlap with base starts at its N term, or the
                         overlap with cand starts at its C term, then we only
                         include cand -> base; and vice versa.
                         */
                        if (baseOverlap > 0 && candOverlap < candChain.residueSize() - numResOverlap) {
                            // Add base -> cand
                            FuseCandidate fuseCandidateRep;
                            fuseCandidateRep.overlapSize = numResOverlap;
                            fuseCandidateRep.setStructure1(seed, baseChain.getID());
                            fuseCandidateRep.overlapPosition1 = baseOverlap;
                            fuseCandidateRep.setStructure2(fuseCandidate, candChain.getID());
                            fuseCandidateRep.overlapPosition2 = candOverlap;
                            fuseCandidateRep.rmsd = get<2>(bestOverlap);
                            result.push_back(fuseCandidateRep);
                        }
                        if (candOverlap > 0 && baseOverlap < baseChain.residueSize() - numResOverlap) {
                            // Add cand -> base
                            FuseCandidate fuseCandidateRep;
                            fuseCandidateRep.overlapSize = numResOverlap;
                            fuseCandidateRep.setStructure1(fuseCandidate, candChain.getID());
                            fuseCandidateRep.overlapPosition1 = candOverlap;
                            fuseCandidateRep.setStructure2(seed, baseChain.getID());
                            fuseCandidateRep.overlapPosition2 = baseOverlap;
                            fuseCandidateRep.rmsd = get<2>(bestOverlap);
                            result.push_back(fuseCandidateRep);
                        }
                        // break;
                    }
                    
                } else {
                    vector<Atom *> candAtoms = candChain.getAtoms();
                    vector<Atom *> baseAtoms = baseChain.getAtoms();
                    
                    FuseCandidate fuseCandidateRep;
                    fuseCandidateRep.overlapSize = numResOverlap;
                    
                    if (searchMode == fromCTerm) {
                        // Delete from N terminus of base
                        for (int n = 0; n < baseChain.residueSize() - numResOverlap; n++) {
                            baseAtoms.erase(baseAtoms.begin(), baseAtoms.begin() + baseChain.getResidue(n).atomSize());
                        }
                        // Delete from C terminus of candidate
                        for (int m = candChain.residueSize() - 1; m >= numResOverlap; m--) {
                            candAtoms.erase(candAtoms.end() - candChain.getResidue(m).atomSize(), candAtoms.end());
                        }
                        fuseCandidateRep.setStructure1(seed, baseChain.getID());
                        fuseCandidateRep.setStructure2(fuseCandidate, candChain.getID());
                        fuseCandidateRep.overlapPosition1 = baseChain.residueSize() - numResOverlap;
                        fuseCandidateRep.overlapPosition2 = 0;
                    } else if (searchMode == fromNTerm) {
                        // Delete from C terminus of base
                        for (int m = baseChain.residueSize() - 1; m >= numResOverlap; m--) {
                            baseAtoms.erase(baseAtoms.end() - baseChain.getResidue(m).atomSize(), baseAtoms.end());
                        }
                        // Delete from N terminus of candidate
                        for (int n = 0; n < candChain.residueSize() - numResOverlap; n++) {
                            candAtoms.erase(candAtoms.begin(), candAtoms.begin() + candChain.getResidue(n).atomSize());
                        }
                        fuseCandidateRep.setStructure1(fuseCandidate, candChain.getID());
                        fuseCandidateRep.setStructure2(seed, baseChain.getID());
                        fuseCandidateRep.overlapPosition1 = candChain.residueSize() - numResOverlap;
                        fuseCandidateRep.overlapPosition2 = 0;
                    }
                    
                    float rmsd = maxDeviation(candAtoms, baseAtoms); // rms.rmsd(candAtoms, baseAtoms);
                    if (rmsd <= rmsdCutoff) {
                        //cout << "RMSD: " << rmsd << endl;
                        fuseCandidateRep.rmsd = rmsd;
                        result.push_back(fuseCandidateRep);
                        break;
                    }
                }
            }
        }
    }
    
    return result;

}

vector<FuseCandidate> FuseCandidateFinder::findOverlappingSeeds(Structure *seed, vector<Structure *> &structures, ProximitySearch &ps) {
    // Find the set of structures that could possibly have an rmsd less than the threshold
    unordered_set<Structure *> filteredCandidates;
    mstreal singleAtomThreshold = sqrt(rmsdCutoff * rmsdCutoff * numResOverlap * 4);
    for (Atom *a: seed->getAtoms()) {
        CartesianPoint coor = a->getCoor();
        vector<int> nearbyTags = ps.getPointsWithin(coor, 0.0, singleAtomThreshold, true);
        for (int tag: nearbyTags) {
            filteredCandidates.insert(structures[tag]);
        }
    }

    vector<Structure *> candVector(filteredCandidates.begin(), filteredCandidates.end());
    return findOverlappingSeeds(seed, candVector);
}

bool FuseCandidateFinder::checkOverlapAlignment(const vector<Atom *> &atoms1, const vector<Atom *> &atoms2) {
    vector<Atom *> alphas1;
    for (Atom *atom: atoms1) {
        if (atom->isNamed("CA"))
            alphas1.push_back(atom);
    }
    vector<Atom *> alphas2;
    for (Atom *atom: atoms2) {
        if (atom->isNamed("CA"))
            alphas2.push_back(atom);
    }
    
    if (alphas1.size() < 2 || alphas2.size() < 2)
        return true;
    
    mstreal dx1 = alphas1[alphas1.size() - 1]->getX() - alphas1[0]->getX();
    mstreal dy1 = alphas1[alphas1.size() - 1]->getY() - alphas1[0]->getY();
    mstreal dz1 = alphas1[alphas1.size() - 1]->getZ() - alphas1[0]->getZ();
    mstreal dx2 = alphas2[alphas2.size() - 1]->getX() - alphas2[0]->getX();
    mstreal dy2 = alphas2[alphas2.size() - 1]->getY() - alphas2[0]->getY();
    mstreal dz2 = alphas2[alphas2.size() - 1]->getZ() - alphas2[0]->getZ();
    mstreal dot = dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
    mstreal mag1 = sqrt(pow(dx1, 2.0) + pow(dy1, 2.0) + pow(dz1, 2.0));
    mstreal mag2 = sqrt(pow(dx2, 2.0) + pow(dy2, 2.0) + pow(dz2, 2.0));
    return dot / (mag1 * mag2) > minCosAngle;
}

mstreal FuseCandidateFinder::maxDeviation(const vector<Atom *> &s1, const vector<Atom *> &s2) {
    MstUtils::assert(s1.size() == s2.size(), "Can't compute max deviation of unequal-sized atom vectors");
    mstreal result = 0.0;
    for (int i = 0; i < s1.size(); i++) {
        Atom *a1 = s1[i];
        Atom *a2 = s2[i];
        MstUtils::assert(a1->getName() == a2->getName(), "Atoms are mismatched");
        mstreal dx = a1->getX() - a2->getX();
        mstreal dy = a1->getY() - a2->getY();
        mstreal dz = a1->getZ() - a2->getZ();
        result = max(result, sqrt(pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0)));
    } 
    return result;
}

mstreal FuseCandidateFinder::atomsFormOverlap(const vector<Atom *> &atoms1, const vector<Atom *> &atoms2, mstreal threshold) {
    if (atoms1.size() != numResOverlap * 4 || atoms2.size() != numResOverlap * 4)
        return DBL_MAX;

    mstreal realThreshold = min(threshold, (mstreal)rmsdCutoff);
    mstreal rmsd = quickRMSD(atoms1, atoms2, realThreshold);

    return (rmsd <= realThreshold && checkOverlapAlignment(atoms1, atoms2)) ? rmsd : DBL_MAX;
}

vector<tuple<int, int, float>> FuseCandidateFinder::scanOverlaps(Chain &baseChain, Chain &candChain) {
    high_resolution_clock::time_point startTime = high_resolution_clock::now();
    RMSDCalculator rms;
    vector<tuple<int, int, float>> results;
    for (int baseStart = 0; baseStart < baseChain.residueSize() - numResOverlap; baseStart++) {
        tuple<int, int, float> bestResult = make_tuple(0, 0, rmsdCutoff);
        for (int candStart = 0; candStart < candChain.residueSize() - numResOverlap; candStart++) {
            // Test these overlaps
            vector<Atom *> candAtoms = candChain.getAtoms();
            vector<Atom *> baseAtoms = baseChain.getAtoms();
            
            // Delete around base
            int deleteIndex = 0;
            for (int n = 0; n < baseChain.residueSize(); n++) {
                if (n >= baseStart && n < baseStart + numResOverlap)
                    deleteIndex += baseChain.getResidue(n).atomSize();
                else
                    baseAtoms.erase(baseAtoms.begin() + deleteIndex, baseAtoms.begin() + deleteIndex + baseChain.getResidue(n).atomSize());
            }
            // Delete around candidate
            deleteIndex = 0;
            for (int n = 0; n < candChain.residueSize(); n++) {
                if (n >= candStart && n < candStart + numResOverlap)
                    deleteIndex += candChain.getResidue(n).atomSize();
                else
                    candAtoms.erase(candAtoms.begin() + deleteIndex, candAtoms.begin() + deleteIndex + candChain.getResidue(n).atomSize());
            }

            double threshold = get<2>(bestResult);
            double rmsd = atomsFormOverlap(candAtoms, baseAtoms, threshold);
            if (rmsd < DBL_MAX) {
                bestResult = make_tuple(baseStart, candStart, rmsd);
            }
        }
        if (get<2>(bestResult) <= rmsdCutoff)
            results.push_back(bestResult);
    }
    high_resolution_clock::time_point endTime = high_resolution_clock::now();
    timeOld += duration_cast<nanoseconds>( endTime - startTime).count();
    timeTrials++;
    return results;
}

void FuseCandidateFinder::writeFuseCandidates(vector<Structure *> seeds, vector<Structure *> candidates, string outDir, string outFileName, string baseSeedPath, bool shouldAppend, int candidatesBatchIdx) {
    FuseCandidateFile outFile(MstSystemExtension::join(outDir, outFileName), shouldAppend);
    
    timeTrials = 0;
    timeNew = 0.0;
    timeOld = 0.0;
    
    // New 3/22/20: Build a proximity search of candidates instead of iterating by brute force
    // The proximity search object will store atoms tagged by the index of the candidate in the candidates list
    vector<mstreal> bbox;
    if (candidatesBatchIdx < 0 || batchBoundingBoxes.count(candidatesBatchIdx) == 0)
        bbox = computeBatchBoundingBox(candidates);
    else
        bbox = batchBoundingBoxes[candidatesBatchIdx];

    ProximitySearch ps(bbox[0], bbox[2], bbox[4], bbox[1], bbox[3], bbox[5], 30);
    for (int candIdx = 0; candIdx < candidates.size(); candIdx++) {
        Structure *c = candidates[candIdx];
        AtomPointerVector apv(c->getAtoms());
        vector<int> indexes(apv.size(), candIdx);
        ps.addAtoms(apv, &indexes);
    }

    for (int j = 0; j < seeds.size(); j++) {
        Structure *seed = seeds[j];
        if (j % 100 == 0)
            cout << j << ": " << MstSystemExtension::fileName(seed->getName()) << endl;

        vector<FuseCandidate> results = findOverlappingSeeds(seed, candidates); //, ps); // );
        outFile.write(results, baseSeedPath);
    }
    cout << "New: " << timeNew / timeTrials << " vs old: " << timeOld / timeTrials << endl;
}

void FuseCandidateFinder::writeFuseCandidates(vector<string> candidatePaths, vector<string> *chainIDs, string outDir, string baseSeedPath) {
    PairStructureIterator structureIter(candidatePaths, 5000, chainIDs);
    int i = 0;
    while (structureIter.hasNext()) {
        if (i % batchModulus == batchModValue) {
            auto batches = structureIter.next();
            cout << "Batch " << i << endl;
            vector<Structure *> seeds = get<0>(batches);
            vector<Structure *> candidates = get<1>(batches);
            writeFuseCandidates(seeds, candidates, outDir, "overlaps" + to_string(batchModValue) + ".txt", baseSeedPath);
        } else {
            structureIter.skip();
        }
        i++;
    }
}

void FuseCandidateFinder::writeFuseCandidates(string seedPath, string seedChainIDs, vector<string> candidatePaths, vector<string> *chainIDs, string outPath, string baseSeedPath) {
    
    Structure seed(seedPath);
    Structure trimmedSeed;
    if (seedChainIDs.empty())
        trimmedSeed = seed;
    else
        extractChains(seed, seedChainIDs, trimmedSeed);
    
    StructureIterator iter(candidatePaths, 10000, chainIDs);
    FuseCandidateFile outFile(outPath);
    
    while (iter.hasNext()) {
        vector<Structure *> candidates = iter.next();
        vector<FuseCandidate> results = findOverlappingSeeds(&seed, candidates);
        outFile.write(results, baseSeedPath);
    }
}

void FuseCandidateFinder::writeFuseCandidates(string binaryFilePath, string outDir) {
    BatchPairStructureIterator structureIter(binaryFilePath, batchModValue, batchModulus, 1000, "0");
    int i = 0;
    bool firstBatch = true;
    batchBoundingBoxes.clear();
    while (structureIter.hasNext()) {
        auto batches = structureIter.next();
        cout << "Batch " << i << endl;
        i++;
        int firstIndex = structureIter.getFirstIndex();
        int secondIndex = structureIter.getSecondIndex();
        vector<Structure *> seeds = batches.first;
        vector<Structure *> candidates = batches.second;

        // Check if batches overlap first
        if (batchBoundingBoxes.count(firstIndex) == 0)
            batchBoundingBoxes[firstIndex] = computeBatchBoundingBox(seeds);
        if (batchBoundingBoxes.count(secondIndex) == 0)
            batchBoundingBoxes[secondIndex] = computeBatchBoundingBox(candidates);

        vector<mstreal> b1 = batchBoundingBoxes[firstIndex];
        vector<mstreal> b2 = batchBoundingBoxes[secondIndex];
        if (!intervalsOverlap(b1[0], b1[1], b2[0], b2[1]) ||
            !intervalsOverlap(b1[2], b1[3], b2[2], b2[3]) ||
            !intervalsOverlap(b1[4], b1[5], b2[4], b2[5])) {
            cout << "Eliminating because bounding boxes are disjoint" << endl;
            continue;
        }

        writeFuseCandidates(seeds, candidates, outDir, "overlaps" + to_string(batchModValue) + ".txt", "", !firstBatch); // Append to the existing file for non-first batches
        firstBatch = false;
    }
    cout << "Done" << endl;
}

/**
 Assumes lo1 < hi1 and lo2 < hi2
 */
bool FuseCandidateFinder::intervalsOverlap(mstreal lo1, mstreal hi1, mstreal lo2, mstreal hi2) {
    return !(hi1 < lo2 || hi2 < lo1);
}
