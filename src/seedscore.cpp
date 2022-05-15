#include <stdio.h>
#include "seedscore.h"
#include "utilities.h"
#include <chrono>

using namespace std::chrono;

string writePerResidueScores(unordered_map<Residue*, mstreal> scores) {
    ostringstream field;
    for (auto it: scores) {
        field << it.first->getName() << it.first->getResidueIndex() << "=" << it.second << ";";
    }
    return field.str();
}

void detailedPrintFragment(Fragment &frag) {
    cout << "{" << endl;
    cout << "\tResidues: ";
    for (Residue *res: frag.getResidues()) {
        cout << res->getChainID() << res->getResidueIndexInChain() << res->getName() << ",";
    }
    cout << endl;
    cout << "\tContacts: ";
    for (Residue *res: frag.getContacts()) {
        cout << res->getChainID() << res->getResidueIndexInChain() << res->getName() << ",";
    }
    cout << endl;
    cout << "\tExpansions:" << endl;
    auto expansions = frag.getExpansions();
    for (Residue *res: frag.getResidues()) {
        cout << "\t\t" << res->getChainID() << res->getResidueIndexInChain() << res->getName() << ": ";
        auto it = expansions.find(res);
        if (it != expansions.end()) {
            for (Residue *expRes: it->second) {
                cout << expRes->getChainID() << expRes->getResidueIndexInChain() << expRes->getName() << ",";
            }
        } else {
            cout << "<none>";
        }
        cout << endl;
    }
    cout << "}" << endl;
}

#pragma mark - FASSTScorer

FASSTScorer::FASSTScorer(Structure *target, string configFilePath, double fractionIdentity, int maxNumMatches, double vdwRadius): SeedScorer(target), fractionIdentity(fractionIdentity), maxNumMatches(maxNumMatches), config(configFilePath), vdwRadius(vdwRadius) {
    loadFromTarget(config.getDB());
}

FASSTScorer::~FASSTScorer() {
    delete fasst;
    delete rl;
}

void FASSTScorer::loadFromTarget(string fasstDB) {
    unordered_map<string, Structure> parentStructs;
    
    // Set up the target structure by removing sidechains
    targetStructBB = Structure(*_target);
    removeSideChains(targetStructBB);
    
    // all this is backbone only
    vector<Residue*> targetResidues = targetStructBB.getResidues();
    targetResidueSet.insert(targetResidues.begin(), targetResidues.end());
    psTargetAPV = targetStructBB.getAtoms();
    ps = ProximitySearch(psTargetAPV, 40.0);
    
    rl = new RotamerLibrary(config.getRL());
    fasst = new FASST();
    fasst->setRedundancyProperty("sim");
    fasst->setMaxNumMatches(maxNumMatches);
    
    if (fasst->numTargets() == 0) {
        cout << "Reading FASST pdbs" << endl;
        fasst->readDatabase(fasstDB, 2); // memory-save flag
        cout << "Done reading FASST pdbs, has " << fasst->numTargets() << " member structures" << endl;
    }
    //cout << "Memory usage: " << MstSys::memUsage() << endl;
}

bool FASSTScorer::prepareCombinedStructure(Structure *seed) {
    // Clean up the seed
    Structure tmpStructure(*seed);
    Structure pose;
    cleanStructure(tmpStructure, pose, true, true, true);
    
    if (pose.residueSize() == 0) {
        cout << "seed has no residues" << endl;
        return false;
    } else if (pose.chainSize() > 1) {
        cout << "seed has more than one chain" << endl;
        return false;
    }
    
    removeSideChains(pose); // may not always want to do this
    AtomPointerVector poseAPV = pose.getAtoms();
    set<int> structResis; // empty since we are not excluding any chains
    int adjustNum = 0;
    
    //cout << "Checking for clash..." << endl;
    if ((!ignoreClash) & (isClash(ps, psTargetAPV, poseAPV, vdwRadius))) {
        cout << "pose clashes with target" << endl;
        return false;
    }
    //cout << "Checked clash" << endl;
    
    // Add the seed to the target as a new chain
    Chain* c = new Chain(generateChainID(targetStructBB.chainSize()), "", pose.getResidues());
    for (int j = 0; j < c->residueSize(); j++) {
        Residue* res = &c->getResidue(j);
        res->setNum(j + 1);
    }
    targetStructBB.appendChain(c); // add chain c to the target
    
    return true;
}

void FASSTScorer::resetCombinedStructure() {
    Chain* poseChain = &targetStructBB[targetStructBB.chainSize() - 1];
    targetStructBB.deleteChain(poseChain);
}

unordered_map<Residue *, mstreal> FASSTScorer::invalidScoreMap(Structure *seed) {
    unordered_map<Residue *, mstreal> result;
    for (Residue *res: seed->getResidues()) {
        result[res] = DBL_MAX;
    }
    return result;
}

unordered_map<Residue*, mstreal> FASSTScorer::remapResiduesFromCombinedStructure(unordered_map<Residue *, mstreal> tempResult, Structure *seed) {
    unordered_map<Residue *, mstreal> result;
    for (auto item: tempResult) {
        result[&seed->getResidue(item.first->getResidueIndexInChain())] = item.second;
    }
    return result;
}

#pragma mark - contactCounter

int contactCounter::countContacts(Structure* seed) {
    int totalNewContacts = 0;
    
    // Add copy of seed chain to target_copy
    if (seed->chainSize() != 1) MstUtils::error("Seed should have a single chain","contactCounter::countContacts");
    Chain* seedChain = new Chain(seed->getChain(0));
    targetCopy->appendChain(seedChain,false);
    
    // Get seed residues, excluding those within flankingRes from the termini
//    vector<Residue*> seedRes = getSeedRes(seedChain);
    vector<Residue*> seedRes = seedChain->getResidues();
    vector<Residue*> seedResToCount = eliminateFlankingRes(seedRes);
    
    if (seedRes.empty()) return totalNewContacts;
    
    if (verbose) cout << "Seed chain with " << seedRes.size() << " residues total, " << seedResToCount.size() << " that will be searched for contacts to the target" << endl;
    
    // Get all kinds of contacts
    contactList bbConts; contactList sbConts; contactList bsConts; contactList ssConts;
    bool intra = false, verbose = false;
    // note from craig: switched bsConts and sbConts here as we want sb to be SC from target
    // normal order bbConts, bsConts, sbConts, ssConts
    splitContacts(*targetCopy, seedResToCount, RL, cParams, intra, bbConts, sbConts, bsConts, ssConts, verbose);
    bool order = false;
    contactList conts = contactListUnion({bbConts, sbConts, bsConts, ssConts}, order);
    
    totalNewContacts = conts.size();
    if (verbose) cout << "Contacts between seed residues and the target: " << totalNewContacts << endl;
    
    map<Residue*,int> seedResContactCounts;
    if (contact_out != nullptr) for (Residue* R : seedRes) seedResContactCounts[R] = 0;
        
    for (int i = 0; i < conts.size(); i++) {
        Residue *sourceResInCopy = conts.srcResidue(i);
        Residue *dstResInCopy = conts.dstResidue(i);
        int idxInChain = sourceResInCopy->getResidueIndexInChain();
        
        if (contact_out != nullptr) {
           // Add to seed residue contact counts
            seedResContactCounts[sourceResInCopy] += 1;
        }
        
        // Find the destination residue in the target and add to contact count
        Residue *dstResInTarget = getTargetResidue(dstResInCopy);
        contactCounts[dstResInTarget] += 1;
    }
    
    if (contact_out != nullptr) {
        for (auto it : seedResContactCounts) {
            *contact_out << seed->getName() << ",";
            *contact_out << it.first->getChainID() << ",";
            *contact_out << it.first->getNum() << ",";
            *contact_out << it.first->getResidueIndexInChain() << ",";
            *contact_out << it.second << endl;
        }
    }
    
    // Reset targetCopy
    targetCopy->deleteChain(seedChain);
    
    return totalNewContacts;
}

void contactCounter::readContactsFile(string contactsFile) {
    vector<string> lines = MstUtils::fileToArray(contactsFile);
    
    for (string line : lines) {
        if (line == lines.front()) continue; // skip the header line
        vector<string> line_split = MstUtils::split(line,"\t");
        if (line_split.size() != 4) MstUtils::error("Line in contactsFile should have 4 fields, but instead has: "+MstUtils::toString(line_split.size()),"contactCounter::readContactsFile");
        string chainID = line_split[0];
        int resNum = MstUtils::toInt(line_split[1]);
        string resName = line_split[2];
        int count = MstUtils::toInt(line_split[3]);
        
        // Check if residue exists in the target protein
        Chain *targetChain = target->getChainByID(chainID);
        if (targetChain == NULL) MstUtils::error("Chain with ID: "+chainID+" not found in target","contactCounter::readContactsFile");
        Residue* targetRes = targetChain->findResidue(resName, resNum);
        if (targetRes == NULL) MstUtils::error("Residue with num,name "+MstUtils::toString(resNum)+","+resName+" not found in target","contactCounter::readContactsFile");
        
        // Add to count
        contactCounts[targetRes] += count;
    }
}

void contactCounter::writeContactsFile(string contactsFile) {
    fstream contacts;
    MstUtils::openFile(contacts,contactsFile,fstream::out);
    string sep = ",";
    
    for (auto it : contactCounts) {
        // get the target residue and necessary info
        Residue *targetRes = it.first;
        string targetResChainID = targetRes->getChainID();
        int resNum = targetRes->getNum();
        string resName = targetRes->getName();
        int count = it.second;
        
        // write to file
        contacts << targetResChainID << sep;
        contacts << resNum << sep;
        contacts << resName << sep;
        contacts << count << endl;
    }
    contacts.close();
}

#pragma mark - SequenceStructureCompatibilityScorer

SequenceStructureCompatibilityScorer::SequenceStructureCompatibilityScorer(Structure *target, rmsdParams& rParams, contactParams& contParams, string configFilePath, int targetFlank, int seedFlank, double fractionIdentity, double minRatio, double pseudocount, int minNumMatches, int maxNumMatches, double vdwRadius): FASSTScorer(target, configFilePath, fractionIdentity, maxNumMatches,   vdwRadius), targetFlank(targetFlank), seedFlank(seedFlank), contParams(contParams), rParams(rParams), minRatio(minRatio), pseudocount(pseudocount), minNumMatches(minNumMatches) {
    fragParams = FragmentParams(max(targetFlank, seedFlank), false);
}

SequenceStructureCompatibilityScorer::~SequenceStructureCompatibilityScorer() {
    if (scoreWriteOut != nullptr)
        delete scoreWriteOut;
}

unordered_map<Residue*, mstreal> SequenceStructureCompatibilityScorer::score(Structure *seed) {
    
    if (!prepareCombinedStructure(seed)) {
        return invalidScoreMap(seed);
    }
    
    Chain* poseChain = &targetStructBB[targetStructBB.chainSize() - 1];
    vector<Residue*> poseResidues = poseChain->getResidues();
    
    // Get contacts
    contactList bbConts; contactList sbConts; contactList bsConts; contactList ssConts;
    bool intra = false;
    // note from craig: switched bsConts and sbConts here as we want sb to be SC from target
    // normal order bbConts, bsConts, sbConts, ssConts
    splitContacts(targetStructBB, poseResidues, rl, contParams, intra, bbConts, sbConts, bsConts, ssConts);
    contactList conts = contactListUnion({bbConts, sbConts, bsConts, ssConts});
    
    if (mustContact && conts.size() == 0) {
        // doesn't contact the correct regions of the target
        cout << "The provided seed doesn't contact the target" << endl;
        resetCombinedStructure();
        return invalidScoreMap(seed);
    }
    
    contactList seqConts = contactListUnion({sbConts, ssConts});
    cout << "Contacts involved with target side chains:" << endl;
    writeContactList(cout, seqConts);
    
    // Perform the scoring
    unordered_map<Residue *, mstreal> combinedResult = sequenceStructureScore(seed, targetStructBB, seqConts, targetResidueSet);
    unordered_map<Residue *, mstreal> result;
    
    if (!countingContacts) {
        // We need to remap the residues from the merged structure to the original seed
        result = remapResiduesFromCombinedStructure(combinedResult, seed);
        cout << "Original result had " << combinedResult.size() << ", new result has " << result.size() << endl;
        
        // Set values in adjacency graph map, if available
        if (adjacencyGraph != nullptr) {
            for (auto pair: result) {
                adjacencyGraph->setValue(pair.first, pair.second);
            }
        }
    }
    
    // Clean up target structure
    resetCombinedStructure();
    
    return result;
}

ofstream *SequenceStructureCompatibilityScorer::getScoreWriteStream() {
    if (scoreWriteOut != nullptr)
        return scoreWriteOut;
    if (scoresWritePath != nullptr) {
        scoreWriteOut = new ofstream(*scoresWritePath);
        if (!scoreWriteOut->is_open()) {
            cerr << "Couldn't open scores write path" << endl;
            scoresWritePath = nullptr;
            scoreWriteOut = nullptr;
            return nullptr;
        }
        
        // Write the header file
        ofstream& out = *scoreWriteOut;
        out << "seed,seed_res_num,target_res_num,num_results,pair_score,bg_score,seq_score,contact_score,total_score" << endl;
    }
    return scoreWriteOut;
}

void SequenceStructureCompatibilityScorer::collectContacts(Structure *seed) {
    countingContacts = true;
    auto result = score(seed);
    countingContacts = false;
}

void SequenceStructureCompatibilityScorer::writeContactCounts(string filePath) {
    ofstream outputSS(filePath, ios::out);
    if (!outputSS.is_open()) {
        cerr << "couldn't open out stream" << endl;
        return;
    }

    for (pair<Residue *, int> item: targetContactCounts) {
        outputSS << item.first->getResidueIndex() << "," << item.second << endl;
    }
}

void SequenceStructureCompatibilityScorer::readContactCounts(string filePath) {
    ifstream readstream(filePath);
    if (!readstream.is_open()) {
        cerr << "couldn't open contact counts file" << endl;
        return;
    }
    
    string line;
    while (getline(readstream, line)) {
        vector<string> comps = splitString(line, ",");
        Residue *targetRes = &targetStructBB.getResidue(atoi(comps[0].c_str()));
        targetContactCounts[targetRes] = atoi(comps[1].c_str());
    }
}

unordered_map<Residue *, mstreal> SequenceStructureCompatibilityScorer::sequenceStructureScore(Structure *seed, Structure &combStruct, contactList &cl, set<Residue *> &toScore) {
    // Store a map from each seed residue to its score. In this method we score
    // residues in the TARGET, but we will add the scores to the appropriate SEED
    // residue.
    unordered_map<Residue*, mstreal> result;

    FragmentParams newFragParams = fragParams;
    newFragParams.pair = true;
    newFragParams.partialFlank = true;
    Fragmenter fragmenter(combStruct, newFragParams);
    fragmenter.fragment(cl);
    set<Fragment> frags = fragmenter.getFragments();
    
    // For each residue to be scored, store a list of fragments around the residue in fragMapping, and the contacting residues themselves in contMap
    map<Residue*, vector<Fragment>> fragMapping;
    map<Residue*, vector<Residue*>> contMap; // residue to be scored in target -> fragment
    for (auto it = frags.begin(); it != frags.end(); it++) {
        Fragment frag = *it;
        vector<Residue*> contact = frag.getID();
        int found1 = toScore.count(contact[0]);
        int found2 = toScore.count(contact[1]);
        
        // In structgen/src/score.cpp, both residues are allowed to be scored. For this use case, I think we only want contacts with one scoreable residue.
        if (found1 + found2 >= 2) {
            cout << "Both residues in contact " + residueID(*contact[0]) + " - " + residueID(*contact[1]) + " cannot be scored. Only one can be scored" << endl;
            continue;
        }
        if (found1 == 1) {
            fragMapping[contact[0]].push_back(frag);
            contMap[contact[0]].push_back(contact[1]);
        } else if (found2 == 1) {
            fragMapping[contact[1]].push_back(frag);
            contMap[contact[1]].push_back(contact[0]);
        }
    }
    cout << "Collected fragments" << endl;
    
    double score = 0;
    map<Fragment, double> seenProbs;
    unordered_set<Residue *> completedSeedResidues; // seed residues with already-cached results from adjacency graph
    ofstream *out = getScoreWriteStream();

    cout << "Iterating over residues" << endl;
    // Iterate over each residue to be scored
    for (auto it = fragMapping.begin(); it != fragMapping.end(); it++) {
        Residue* res = it->first;
        vector<Fragment> frags = it->second;
        vector<Residue*>& resContacts = contMap[res];
        string resCode = SeqTools::tripleToSingle(res->getName());
        cout << "Residue code is " << resCode << endl;
        int resType = SeqTools::seqToIdx(resCode, "")[0];
        // note from craig: should probably use pseudocounts here instead?
        if (resType == SeqTools::unknownIdx()) {
            cout << "Residue type " << res->getName() << " is non standard." << endl;
            for (int i = 0; i < frags.size(); i++) {
                cout << "\t" << "Fragment: ";
                detailedPrintFragment(frags[i]);
            }
            continue;
        }
        
        if (countingContacts) {
            // Increment the number of contacts observed for this target residue
            targetContactCounts[&targetStructBB.getResidue(res->getResidueIndexInChain())] += frags.size();
        } else {
            // Iterate over contacts with this scoring residue, and their associated fragments
        for (int i = 0; i < frags.size(); i++) {
            Fragment& frag = frags[i];
            Residue *seedResidue = resContacts[i];
            if (result[seedResidue] == DBL_MAX || completedSeedResidues.count(seedResidue) != 0)
                continue;
            if (adjacencyGraph != nullptr) {
                mstreal *val = adjacencyGraph->value(&seed->getResidue(seedResidue->getResidueIndexInChain()));
                if (val != nullptr) {
                    cout << "Using cached value for seed residue score" << endl;
                    completedSeedResidues.insert(seedResidue);
                    result[seedResidue] = *val;
                    continue;
                }
            }
            
            if (out != nullptr)
                *out << seed->getName() << "," << seedResidue->getChainID() << seedResidue->getNum() << "," <<  res->getChainID() << res->getNum() << ",";
            
            cout << "Calculating score component" << endl;
            double score = sequenceStructureScoreComponent(frag, res, seedResidue, &seenProbs);
            seenProbs[frag] = score;
            
            cout << "Finished calculating score component" << endl;
            // Add this score to the contacting residue on the SEED chain
            if (score == DBL_MAX) {
                result[seedResidue] = score;
                if (out != nullptr)
                    *out << score << "," << DBL_MAX << "," << DBL_MAX << endl;
            } else {
                result[seedResidue] += score;
                
                // Also add the contact score
                double ctScore = contactScore(seedResidue, res);
                cout << "Scores: " << score << "," << ctScore << endl;
                if (out != nullptr)
                    *out << score << "," << ctScore << "," << score + ctScore << endl;
                result[seedResidue] += ctScore;
            }
        }
        }
    }
    
    if (!countingContacts) {
    uniqueResiduesScored += result.size() - completedSeedResidues.size();
    residuesScored += result.size();
    }
    
    return result;
}

vector<Residue *> SequenceStructureCompatibilityScorer::trimFragment(Fragment &frag, Residue *scoringRes, Residue *seedRes, Structure &result) {
    vector<Residue *> residues;
    int scoringResIndex = scoringRes->getResidueIndex();
    int seedResIndex = seedRes->getResidueIndex();
    vector<Residue *> fragResidues = frag.getResidues();
    for (Residue *res: fragResidues) {
        if (abs(res->getResidueIndex() - scoringResIndex) <= targetFlank || abs(res->getResidueIndex() - seedResIndex) <= seedFlank) {
            residues.push_back(res);
        }
    }
    //MstUtils::assert(find(residues.begin(), residues.end(), seedRes) == residues.end(), "isolating fragment already contains seed residue");
    //residues.push_back(seedRes);
    result = Structure(residues);
    result.setName(frag.toString());
    return residues;
}

mstreal SequenceStructureCompatibilityScorer::sequenceStructureScoreComponent(Fragment &frag, Residue *scoringRes, Residue *seedRes, map<Fragment, double> *seenProbs, bool cacheBackground) {
    cout << "Analyzing contact " << *scoringRes << " - " << *seedRes << " with " << frag.numResidues() << " residues in fragment" << endl;
    string resCode = SeqTools::tripleToSingle(scoringRes->getName());
    int resType = SeqTools::seqToIdx(resCode, "")[0];
    ofstream *out = getScoreWriteStream();
    
    // Search the interface-TERM if not seen yet
    double pairProb = 0;
    if (seenProbs != nullptr && seenProbs->count(frag) > 0) {
        pairProb = (*seenProbs)[frag];
        if (out != nullptr)
            *out << INT_MAX << ",";
    } else {
        if (!noQueries) {
            Structure fragStruct;
            vector<Residue*> fragResidues = trimFragment(frag, scoringRes, seedRes, fragStruct);
            cout << "(trimmed fragment has " << fragStruct.residueSize() << " residues)" << endl;
            
            // Write query to path if debugging
            if (queryWritePath != nullptr) {
                string writePath = MstSystemExtension::join(*queryWritePath, "query_" + scoringRes->getName() + to_string(scoringRes->getResidueIndex()) + "_" + seedRes->getName() + to_string(seedRes->getResidueIndex()) + "_target+seed.pdb");
                fragStruct.writePDB(writePath);
            }
            
            fasst->setQuery(fragStruct);
            double rmsdCut = rParams.rmsdCutoff(fragResidues); // could also include combStruct here
            //cout << "\tQuery has " << fragStruct.atomSize() << " atoms, RMSD " << rmsdCut << endl;
            
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            fasst->search();
            //cout << "\tFound " << fasst->numMatches() << " matches below " << rmsdCut << " Angstroms\n";
            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            //cout << "\tSeed search took " << duration_cast<seconds>( endTime - startTime ).count() << " seconds" << endl;
            
            if (out != nullptr)
                *out << fasst->numMatches() << ",";
            
            if (fasst->numMatches() < minNumMatches) {
                cout << "> Below min number of matches " << minNumMatches << endl;
                return DBL_MAX;
            }
            
            int resi = distance(fragResidues.begin(), find(fragResidues.begin(), fragResidues.end(), scoringRes));
            MstUtils::assert(resi < fragResidues.size(), "residue index out of bounds");
            fasstSolutionSet solSet = fasst->getMatches();
            
            // Generates a 20 x L matrix, where L is length of the sequences
            Matrix freqs = seqFreqs(fasst->getMatchSequences(solSet, FASST::matchType::REGION), pseudocount, true);
            pairProb = freqs(resType, resi);
        }
        
        (*seenProbs)[frag] = pairProb;
        seedQueries++;
    }
    
    // Get the flanking region around the scoring residue
    vector<Residue*> segResidues = frag.getExpansion(scoringRes);
    auto backProbIt = backgroundProbs.find(scoringRes);
    double backProb = -1;
    if (backProbIt != backgroundProbs.end() && cacheBackground) {
        // Use cached value
        cout << "Using cached value for background probability" << endl;
        if (backProbIt->second.first != segResidues.size()) {
            cout << "> Cached probability has a different number of flanking residues! You might want to fix this." << endl;
        }
        backProb = backProbIt->second.second;
    } else {
        if (!noQueries) {
            // Search for the flanking region
            Structure segStruct(segResidues);
            
            double backRMSDCut = rParams.rmsdCutoff(segResidues);
            if (queryWritePath != nullptr)
                segStruct.writePDB(MstSystemExtension::join(*queryWritePath, "query_" + scoringRes->getName() + to_string(scoringRes->getResidueIndex()) + "_target.pdb"));
            fasst->setQuery(segStruct);
            high_resolution_clock::time_point startTime = high_resolution_clock::now();
            fasstSolutionSet segSolSet = fasst->search();
            high_resolution_clock::time_point endTime = high_resolution_clock::now();
            //cout << "\tBackground search took " << duration_cast<seconds>( endTime - startTime ).count() << " seconds" << endl;
            
            // Find index of the scoring residue in the flanking region
            auto resiIt = find(segResidues.begin(), segResidues.end(), scoringRes);
            if (resiIt == segResidues.end()) {
                cout << "> Warning: scoring residue not found in its own fragment" << endl;
                return DBL_MAX;
            }
            int segResi = distance(segResidues.begin(), resiIt);
            
            vector<Sequence> segSeqs = fasst->getMatchSequences(segSolSet, FASST::matchType::REGION);
            Matrix segFreqs = seqFreqs(segSeqs, pseudocount, true);
            backProb = segFreqs(resType, segResi);
            int backNumMatches = fasst->numMatches();
            //cout << "\tFound " << backNumMatches << " background matches below " << backRMSDCut << " Angstroms\n";
        } else {
            backProb = 0.001;
        }
        
        if (cacheBackground)
            backgroundProbs[scoringRes] = make_pair(segResidues.size(), backProb);
        
        targetQueries++;
    }
    
    if (out != nullptr)
        *out << pairProb << "," << backProb << ",";

    // Compute the score for this contact
    double ratio = pairProb / backProb;
    if (ratio < minRatio) {
        cout << "> Ratio " << ratio << " is below " << minRatio << endl;
        return DBL_MAX;
    }
    return -log(ratio);
}

mstreal SequenceStructureCompatibilityScorer::contactScore(Residue *seedRes, Residue *targetRes, double pseudocount) {
    if (targetContactCounts.size() == 0 || adjacencyGraph == nullptr)
        return 0.0;
    
    Residue *realTargetRes = &targetStructBB.getResidue(targetRes->getResidueIndexInChain());
    double numerator = (double)adjacencyGraph->neighborhood(seedRes).residueSize() + pseudocount;
    double denominator = (double)targetContactCounts[realTargetRes] + pseudocount;
    
    return -log(numerator / denominator);
}

#pragma mark - StructureCompatibilityScorer

//SeedDesignabilityScorer::SeedDesignabilityScorer(Structure *target, FragmentParams& fragParams, rmsdParams& rParams, contactParams& contParams, string configFilePath, double fractionIdentity, int minNumMatches, int maxNumMatches, double vdwRadius, bool _scoreAll)

unordered_map<Residue*, mstreal> SeedDesignabilityScorer::score(Structure *seed) {
    _numDesignable = 0;

    // Stores combined structure in targetStructBB
    if (!prepareCombinedStructure(seed)) {
        return invalidScoreMap(seed);
    }
    //cout << "Prepared combined structure" << endl;
    
    Chain* poseChain = &targetStructBB[targetStructBB.chainSize() - 1];
    vector<Residue*> poseResidues = poseChain->getResidues();
    //cout << "Got pose residues" << endl;
    
    contactList bbConts; contactList sbConts; contactList bsConts; contactList ssConts;
    // note from craig: switched bsConts and sbConts here as we want sb to be SC from target
    // normal order bbConts, bsConts, sbConts, ssConts
    splitContacts(targetStructBB, poseResidues, rl, contParams, false, bbConts, sbConts, bsConts, ssConts);
    //cout << "Split contacts" << endl;
    contactList conts = contactListUnion({bbConts, sbConts, bsConts, ssConts});
    //cout << "Took contact list union" << endl;
    
    if (mustContact && conts.size() == 0) {
        // doesn't contact the correct regions of the target
        resetCombinedStructure();
        return invalidScoreMap(seed);
    }
    
    //cout << "Contacts:" << endl;
    writeContactList(cout, conts);
    
    // Compute designability score using the list of contacts
    auto tempResult = designabilityScore(targetStructBB, conts, poseResidues);
    auto result = remapResiduesFromCombinedStructure(tempResult, seed);
    
    // Clean up target structure
    resetCombinedStructure();
    
    return result;
}

bool SeedDesignabilityScorer::clashes(Structure *seed) {
    if (!prepareCombinedStructure(seed))
        return true;
    resetCombinedStructure();
    return false;
}

void SeedDesignabilityScorer::score(Structure *seed, mstreal &totalScore, int &numContacts, int &numDesignable, bool intra, bool score_all) {
    _numDesignable = 0;
    numContacts = 0;
    numDesignable = 0;
    totalScore = DBL_MAX;
      
    // Stores combined structure in targetStructBB
    if (!prepareCombinedStructure(seed)) {
        return;
    }
    
    cout << "Prepared combined structure..." << endl;
    
    Chain* poseChain = &targetStructBB[targetStructBB.chainSize() - 1];
    vector<Residue*> poseResidues = poseChain->getResidues();
    //cout << "Got pose residues" << endl;
    
    contactList bbConts; contactList sbConts; contactList bsConts; contactList ssConts;
    // note from craig: switched bsConts and sbConts here as we want sb to be SC from target
    // normal order bbConts, bsConts, sbConts, ssConts
    splitContacts(targetStructBB, poseResidues, rl, contParams, intra, bbConts, sbConts, bsConts, ssConts);
    contactList conts = contactListUnion({bbConts, sbConts, bsConts, ssConts});
    numContacts = conts.size();
    
    totalScore = 0.0;
    
    if (mustContact && conts.size() == 0) {
        cout << "No contacts to score" << endl;
        resetCombinedStructure();
        return;
    }
    
    writeContactList(cout, conts);
  
    // Compute designability score using the list of contacts
    unordered_map<Residue *, mstreal> tempResult = designabilityScore(targetStructBB, conts, poseResidues, seed->getName());
    for (auto item: tempResult) {
        totalScore += item.second;
    }
    /* Added by sebastian on 20/07/20:
     total score should be normalized by number of residues
    */
    totalScore /= poseResidues.size();
    numDesignable = _numDesignable;
    
    // Clean up target structure
    resetCombinedStructure();
    
}

unordered_map<Residue*, mstreal> SeedDesignabilityScorer::designabilityScore(Structure &combStruct, contactList &cl, vector<Residue*> seedResidues, string seedName) {
    ofstream* out = getScoreWriteStream();
    
    // Fragment structure based on contacts
    FragmentParams newFragParams = fragParams;
    newFragParams.pair = true;
    newFragParams.partialFlank = true;
    Fragmenter fragmenter(combStruct, newFragParams);
    fragmenter.fragment(cl);
    set<Fragment> frags = fragmenter.getFragments();
    
    unordered_map<Residue*, mstreal> result;
    for (Residue* res: seedResidues)
        result[res] = 0;
    
    if (frags.size() == 0) {
        // No fragments - it's trivially designable
        return unordered_map<Residue*, mstreal>();
    }
    
    int prevSufficientNumMatches = fasst->getSufficientNumMatches();
    /**
        If the homology cutoff is in effect, we can't stop searching at minNumMatches, because some of the matches are later filtered
     due to homology with the query, then we would not have enough in our solution set for the fragment to be considered designable.
     Theoretically we should find all matches, but in practice since the DB is always pre-filtered by sequence similarity, we expect no more
     than one match to be filtered out. To be safe, we try to find an extra five matches.
     */
    if (queryHomologyCutoff <= 1.0) fasst->setSufficientNumMatches(minNumMatches+5);
    else fasst->setSufficientNumMatches(minNumMatches);
    
    cout << "searching for matches to " << frags.size() << " fragments " << endl;
    MstTimer timer;
    for (Fragment frag: frags) {
        cout << "Fragment with name: " << frag << endl;
        Structure fragStructure = frag.getStructure();
        
        // Get seed residues, which will be scored
        vector<Residue*> residuesToScore;
        vector<Residue*> contactPair = frag.getID();
        copy_if(contactPair.begin(), contactPair.end(), back_inserter(residuesToScore), [&](Residue* res) { return find(seedResidues.begin(), seedResidues.end(), res) != seedResidues.end() && result[res] != DBL_MAX; });
        if (residuesToScore.size() == 0) {
            cout << "Fragment had no seed residues, skipping" << endl;
            continue;
        }
        
        // if the homologyCut value is <= 1.0, will filter out homologous matches
        termData tD;
        if (queryHomologyCutoff <= 1.0) tD.define(contactPair, fragParams.maxNumFlank);
        
        fasst->setQuery(fragStructure);
        double rmsdCut = rParams.rmsdCutoff(fragStructure);
        fasst->setRMSDCutoff(rmsdCut);
        timer.start();
        cout << "Searching, " << fragStructure.atomSize() << " atoms, " << fragStructure.residueSize() << " residues, RMSD " << rmsdCut << endl;
        fasstSolutionSet sols = fasst->search();
        timer.stop();
        
        if (queryHomologyCutoff <= 1.0) {
            cout << "Done searching. Found " << sols.size() << " matches. Took " << timer.getDuration() << "s" << " ... now filtering homologous matches at a cutoff of " << queryHomologyCutoff << endl;
            tD.setMatches(sols, queryHomologyCutoff, fasst);
            cout << "Done filtering. Now have " << tD.numMatches() << endl;
        } else {
            tD.setMatches(sols, queryHomologyCutoff, NULL);
            cout << "Done searching. Found " << tD.numMatches() << " matches. Took " << timer.getDuration() << "s" << endl;
        }
                
        // write to ofstream, if provided
        if (out != nullptr ) {
            vector<Residue*> centralResID = frag.getID();
            if (centralResID.size() != 2) MstUtils::error("Interface fragment has wrong number of center residues","SeedDesignabilityScorer::designabilityScore");
            *out << seedName << ",";
            *out << centralResID[0]->getChainID() << "," << centralResID[0]->getNum() << ",";
            *out << centralResID[1]->getChainID() << "," << centralResID[1]->getNum() << ",";
            *out << rmsdCut << "," << timer.getDuration() << ",";
        }
        if (tD.numMatches() < minNumMatches) {
            if (out != nullptr) *out << 0 << endl;
            if (!scoreAll) {
                // Set the seed residue score to infinity
                for (Residue *res: residuesToScore) {
                    result[res] = DBL_MAX;
                }
                break;
            }
        }
        else {
            if (out != nullptr) *out << 1 << endl;
            _numDesignable++;
            for (Residue *res: residuesToScore) {
                result[res] += 1;
            }
        }
    }
    cout << "took " << timer.getDuration() << "s to find matches to all fragments" << endl;
    fasst->setSufficientNumMatches(prevSufficientNumMatches);
    return result;
}

ofstream *SeedDesignabilityScorer::getScoreWriteStream() {
    if (scoreWriteOut != nullptr)
        return scoreWriteOut;
    if (scoresWritePath != "") {
        scoreWriteOut = new ofstream(scoresWritePath);
        if (!scoreWriteOut->is_open()) {
            cerr << "Couldn't open scores write path" << endl;
            scoresWritePath = nullptr;
            scoreWriteOut = nullptr;
            return nullptr;
        }
        
        // Write the header file
        ofstream& out = *scoreWriteOut;
        out << "seed,seed_res_chain_id,seed_res_num,target_res_chain_id,target_res_num,rmsd,search_time,match" << endl;
    }
    return scoreWriteOut;
}
