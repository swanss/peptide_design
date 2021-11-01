#include "utilities.h"

using namespace MST;

/* --------- configFile --------- */
configFile::configFile(string configFile) {
    vector<string> lines = MstUtils::fileToArray(configFile);
    for (int i = 0; i < lines.size(); i++) {
        string line = MstUtils::trim(MstUtils::removeComment(lines[i], "#"));
        if (line.empty()) continue;
        vector<string> ents = MstUtils::trim(MstUtils::split(line, "="));
        if (ents.size() != 2) MstUtils::error("could not parse parameter line '" + lines[i] + "' from file " + configFile, "configFile::configFile(const string&)");
        if (ents[0].compare("fasstdb") == 0) {
            fasst_DB_path = ents[1];
        } else if (ents[0].compare("rotlib") == 0) {
            RL_path = ents[1];
        } else {
            cout << "unknown parameter name: " << ents[0] << "... ignoring..." << endl;
        }
    }
}


/* --------- generalUtilities --------- */
string generalUtilities::selectionStringFromContactingRes(contactList conts) {
    vector<Residue*> src = conts.srcResidues();
    vector<Residue*> dest = conts.destResidues();
    
    vector<Residue*> src_unique;
    for (Residue* R : src) if (find(src_unique.begin(),src_unique.end(),R) == src_unique.end()) src_unique.push_back(R);
    
    vector<Residue*> dest_unique;
    for (Residue* R : dest) if (find(dest_unique.begin(),dest_unique.end(),R) == dest_unique.end()) dest_unique.push_back(R);
    
    vector<Residue*> all = MstUtils::setunion(src_unique,src_unique);
    return generalUtilities::selectionStringFromRes(all);
}

string generalUtilities::selectionStringFromRes(vector<Residue*> residues) {
    string selection = "";
    for (Residue* res : residues) {
        selection += "(chain " + res->getChainID() + " and resid " + MstUtils::toString(res->getNum()) + ")";
        if (res != residues.back()) selection += " or ";
    }
    return selection;
}

void generalUtilities::copyAtomCoordinates(Structure *S_dest, const Structure *S_source) {
    vector<Atom*> atoms_source = S_source->getAtoms();
    copyAtomCoordinates(S_dest,atoms_source);
}

void generalUtilities::copyAtomCoordinates(Structure *S_dest, vector<Atom*> atoms_source) {
    vector<Atom*> atoms_dest = S_dest->getAtoms();
    if (atoms_dest.size() != atoms_source.size()) {
        string dest_size = MstUtils::toString(atoms_dest.size());
        string source_size = MstUtils::toString(atoms_source.size());
        MstUtils::error("The two objects have a different number of atoms. "+dest_size+" and "+source_size);
    }
    for (int i = 0; i < atoms_dest.size(); i++) {
        atoms_dest[i]->setCoor(atoms_source[i]->getCoor());
    }
}

vector<Structure*> generalUtilities::readStructuresBin(string binFile) {
    vector<Structure*> structures;
    fstream ifs; MstUtils::openFile(ifs, binFile, fstream::in | fstream::binary, "generalUtilities::readStructuresBin");
    while (ifs.peek() != EOF) {
        Structure* S = new Structure();
        S->readData(ifs);
        structures.push_back(S);
    }
    return structures;
}

vector<Structure*> generalUtilities::readStructuresList(string listFile) {
    vector<Structure*> structures;
    vector<string> structure_files = MstUtils::fileToArray(listFile);
    for (int i = 0; i < structure_files.size(); i++) {
        Structure* S = new Structure(structure_files[i]);
        structures.push_back(S);
    }
    return structures;
}

void generalUtilities::writeStructuresBin(vector<Structure*> structures, string binFile) {
    fstream ofs; MstUtils::openFile(ofs, binFile, fstream::out | fstream::binary, "generalUtilities::writeStructuresBin");
    for (int i = 0; i < structures.size(); i++) {
        if (structures[i] == NULL) MstUtils::error("At least one of the structures is not populated","generalUtilities::writeStructuresBin()");
        Structure* S = structures[i];
        S->writeData(ofs);
    }
    ofs.close();
}

void generalUtilities::writeStructuresBin(vector<Structure> structures, string binFile) {
    fstream ofs; MstUtils::openFile(ofs, binFile, fstream::out | fstream::binary, "generalUtilities::writeStructuresBin");
    for (int i = 0; i < structures.size(); i++) {
        Structure& S = structures[i];
        S.writeData(ofs);
    }
    ofs.close();
}

// copied from dTERMen::getContactsWith
// the first position corresponds to the original residue, the second correspond to its contacts
vector<pair<Residue*, Residue*>> generalUtilities::getContactsWith(const vector<Residue*>& source, ConFind& C, int type, mstreal cd_threshold, mstreal int_threshold, mstreal bbInteraction_cutoff, bool verbose) {
    
    set<Residue*> sourceSet = MstUtils::contents(source);
    vector<pair<Residue*, Residue*>> conts;
    set<pair<Residue*, Residue*>> contsSet;
    pair<Residue*, Residue*> c;
    contactList contList;
    if (verbose) {
        cout << "identifying contacts with:";
        for (int i = 0; i < source.size(); i++) cout << " " << *(source[i]);
        cout << endl;
    }
    
    // get all contacts involving the source residues
    for (int cType = 0; cType < 3; cType++) {
        if (cType == 0) {
            contList = C.getContacts(source, cd_threshold);
        }
        else if (cType == 1) {
            contList = C.getInterference(source, int_threshold);
        }
        else if (cType == 2) {
            //don't include res that are already included by default as flank
            contList = C.getBBInteraction(source, bbInteraction_cutoff);
        }
        // go through each and insert into list, in the right order, if it qualifies
        contList.sortByDegree();
        for (int i = 0; i < contList.size(); i++) {
            //check if either of the residues are present in the source set
            bool isInA = (sourceSet.find(contList.residueA(i)) != sourceSet.end());
            bool isInB = (sourceSet.find(contList.residueB(i)) != sourceSet.end());
            
            //based on the selected type, determines whether this contact can be added
            //type 0 - at least one residue should not be in source set
            //type 1 - both must be in source set
            //type 2 - no residue may be in source set (?)
            if (((type == 0) && (isInA == isInB)) || ((type == 1) && !(isInA && isInB)) || ((type == 2) && !(isInA || isInB))) continue;
            
            //the first residue of the pair is from the source set
            if (!isInA) c = pair<Residue*, Residue*>(contList.residueB(i), contList.residueA(i));
            else c = pair<Residue*, Residue*>(contList.residueA(i), contList.residueB(i));
            
            if (verbose) {
                if (cType == 0) cout << "\t" << "contact-degree ";
                else if (cType == 1) cout << "\t" << "interference ";
                else if (cType == 2) cout << "\t" << "backbone-interaction ";
                cout << "contact with " << *(c.second) << " (from " << *(c.first) << "); " << contList.degree(i) << endl;
            }
            
            // add if new contacting pair
            if (contsSet.find(c) == contsSet.end()) {
                contsSet.insert(c);
                conts.push_back(c);
            }
        }
    }
    if (verbose) cout << "Fragmenter::getContactsWith() -> in the end, found " << conts.size() << " contacts" << endl;
    
    return conts;
}

vector<Residue*> generalUtilities::getContactingResidues(vector<pair<Residue*,Residue*>> cont_pairs) {
    /* Get contacts on non-source side including sidechain-sidechain / sidechain backbone / backbone backbone. */
    vector<Residue*> contacting_res;
    for (int i = 0; i < cont_pairs.size(); i++) {
        Residue* R = cont_pairs[i].second;
        if (find(contacting_res.begin(),contacting_res.end(),R) == contacting_res.end()) contacting_res.push_back(R);
    }
    return contacting_res;
}

contactList generalUtilities::getContactsWith(const vector<Residue *>& source, ConFind& C, mstreal threshold, string type) {
    set<Residue*> sourceSet = MstUtils::contents(source);
    contactList contList;
    contactList filtered_contList;
    
    cout << "identifying contacts with:";
    for (int i = 0; i < source.size(); i++) cout << " " << *(source[i]);
    cout << endl;
    
    if (type == "contact") {
        cout << "contact-degree\t" << threshold << endl;
        contList = C.getContacts(source, threshold);
    } else if (type == "interfering") {
        cout << "interfering\t" << threshold << endl;
        contList = C.getInterfering(source, threshold);
    } else if (type == "interfered") {
        cout << "interfered\t" << threshold << endl;
        contList = C.getInterference(source, threshold);
    } else if (type == "bbinteraction") {
        cout << "bbinteraction\t" << threshold << endl;
        contList = C.getBBInteraction(source, threshold); //it's actually a cutoff
    }
    
    // go through each and insert into list, in the right order, if it qualifies
    contList.sortByDegree();
    for (int i = 0; i < contList.size(); i++) {
        //check if either of the residues are present in the source set
        bool A_in_source = (sourceSet.find(contList.residueA(i)) != sourceSet.end());
        bool B_in_source = (sourceSet.find(contList.residueB(i)) != sourceSet.end());
        
        //This assumes A is always the from the source set
        
        //1) A must be in the source set
        //2) B must not be in the source set
        
        if (!(A_in_source) || (B_in_source)) continue;
        
        //the first residue of the pair is from the source set
        cout << "contact with " << *(contList.residueB(i)) << " (from " << *(contList.residueA(i)) << "); " << contList.degree(i) << endl;
        filtered_contList.addContact(contList.residueA(i), contList.residueB(i), contList.degree(i));
        
    }
    cout << "Fragmenter::getContactsWith() -> in the end, found " << filtered_contList.size() << " contacts between source and non-source residues" << endl;
    return filtered_contList;
}

set<pair<Residue*,Residue*>> generalUtilities::mergeContactLists(vector<contactList> CL) {
    set<pair<Residue*,Residue*>> merged_contacts;
    for (int i = 0; i < CL.size(); i++) {
        for (int j = 0; j < CL[i].size(); j++) {
            Residue* R_i = CL[i].residueA(j);
            Residue* R_j = CL[i].residueB(j);
            pair<Residue*,Residue*> cont(R_i,R_j);
            merged_contacts.insert(cont);
        }
    }
    return merged_contacts;
}

// recursive algorithm for generating all unique combinations of k residues from a larger set (without replacement).
vector<vector<Residue*>> generalUtilities::generateAllCombinationsKRes(vector<Residue*> initial_set, int k) {
    vector<vector<Residue*>> final_residue_combinations;
    if (k == 0) {
        vector<vector<Residue*>> remaining_res;
        return remaining_res;
    }
    //decrement k; move down one level of the tree
    k = k - 1;
    
    for (int res_selector = 0; res_selector < initial_set.size(); res_selector++) {
        //extract the "current" residue and remove from initial set
        Residue* R = initial_set[res_selector];
        vector<Residue*> reduced_set(initial_set.begin()+res_selector+1,initial_set.end());
        
        //call function to generate the next layer of the tree
        vector<vector<Residue*>> residue_combinations_suffixes;
        residue_combinations_suffixes = generateAllCombinationsKRes(reduced_set,k);
        
        //for each individual suffix below the current branch, add to "current" residue
        for (int suffix_selector = 0; suffix_selector < residue_combinations_suffixes.size(); suffix_selector++){
            vector<Residue*> new_combination;
            new_combination.reserve(residue_combinations_suffixes[suffix_selector].size() + 1);
            new_combination.insert(new_combination.end(), R);
            new_combination.insert(new_combination.end(), residue_combinations_suffixes[suffix_selector].begin(), residue_combinations_suffixes[suffix_selector].end());
            final_residue_combinations.push_back(new_combination);
        }
        
        //if bottom of tree, then simply add the "current" residue to final residue combinations
        if (( k == 0 ) && (residue_combinations_suffixes.size() == 0)) {
            vector<Residue*> new_combination = {R};
            final_residue_combinations.push_back(new_combination);
        }
    }
    return final_residue_combinations;
};

vector<vector<Residue*>> generalUtilities::generateAllCombinationsRes(vector<Residue*> initial_set) {
    vector<vector<Residue*>> all_combinations;
    for (int k = initial_set.size(); k > 0; k--) {
        vector<vector<Residue*>> all_k_res_combinations = generalUtilities::generateAllCombinationsKRes(initial_set,k);
        all_combinations.insert(all_combinations.end(),all_k_res_combinations.begin(),all_k_res_combinations.end());
    }
    return all_combinations;
}

mstreal generalUtilities::fragDegreesOfFreedom(const vector<int>& L, mstreal L0) {
    mstreal a = (mstreal) exp(-1./L0);
    int N = 0, n;
    mstreal c = 0;
    
    // disjoint segments are counted as independent, so their correlation
    // with respect to each other is zero
    for (int i = 0; i < L.size(); i++) {
        N += L[i];
        n = L[i];
        c = c + (a/(1-a))*(n-1) - pow((a/(1-a)), 2)*(1 - pow(a, n-1));
    }
    double df = N*(1 - (2.0/(N*(N-1)))*c);
    
    //    return rmsdMax/sqrt(N/df);
    return df;
}


mstreal generalUtilities::fragDegreesOfFreedom(const vector<vector<int> >& I, mstreal L0) {
    int N = 0;
    mstreal c = 0;
    
    // disjoint segments are counted as independent, so their correlation
    // with respect to each other is zero
    for (int i = 0; i < I.size(); i++) {
        for (int j = 0; j < I[i].size(); j++) {
            for (int k = j + 1; k < I[i].size(); k++) {
                c = c + exp(-abs(I[i][j] - I[i][k])/L0);
            }
        }
        N += I[i].size();
    }
    mstreal df = N*(1 - (2.0/(N*(N-1)))*c);
    
    //    return rmsdMax/sqrt(N/df);
    return df;
}

mstreal generalUtilities::fragDegreesOfFreedom(const Structure& S, mstreal L0) {
    vector<int> L(S.chainSize());
    for (int i = 0; i < S.chainSize(); i++) L[i] = S[i].residueSize();
    return generalUtilities::fragDegreesOfFreedom(L, L0);
}

mstreal generalUtilities::fragDegreesOfFreedom(const vector<int>& J, const Structure& S, mstreal L0) {
    vector<vector<int> > I;
    map<Chain*, vector<int> > residuesFromChain;
    for (int i = 0; i < J.size(); i++) residuesFromChain[S.getResidue(J[i]).getChain()].push_back(J[i]);
    vector<Chain*> keys = MstUtils::keys(residuesFromChain);
    I.resize(keys.size());
    for (int i = 0; i < keys.size(); i++) I[i] = residuesFromChain[keys[i]];
    return generalUtilities::fragDegreesOfFreedom(I, L0);
}

mstreal generalUtilities::cosAngleBetweenNormalVectors(Residue *R1, Residue *R2) {
    //get the N,Ca,and C backbone atoms for both residues
    vector<Atom*> bb1 = RotamerLibrary::getBackbone(R1);
    vector<Atom*> bb2 = RotamerLibrary::getBackbone(R2);
    if ((bb1[RotamerLibrary::bbCA] == NULL) || (bb1[RotamerLibrary::bbC] == NULL) || (bb1[RotamerLibrary::bbN] == NULL)) {
        MstUtils::error("cannot place rotamer in residue " + MstUtils::toString(*R1) + ", as it lacks proper backbone", "RotamerLibrary::placeRotamer");
    }
    if ((bb2[RotamerLibrary::bbCA] == NULL) || (bb2[RotamerLibrary::bbC] == NULL) || (bb2[RotamerLibrary::bbN] == NULL)) {
        MstUtils::error("cannot place rotamer in residue " + MstUtils::toString(*R2) + ", as it lacks proper backbone", "RotamerLibrary::placeRotamer");
    }
    CartesianPoint CA1 = CartesianPoint(bb1[RotamerLibrary::bbCA]);
    CartesianPoint C1 = CartesianPoint(bb1[RotamerLibrary::bbC]);
    CartesianPoint N1 = CartesianPoint(bb1[RotamerLibrary::bbN]);
    CartesianPoint CA2 = CartesianPoint(bb2[RotamerLibrary::bbCA]);
    CartesianPoint C2 = CartesianPoint(bb2[RotamerLibrary::bbC]);
    CartesianPoint N2 = CartesianPoint(bb2[RotamerLibrary::bbN]);
    
    //take the cross product (N -> Ca) x (Ca -> C) to find the normal vector to the plane formed by N,Ca,C
    CartesianPoint N1toCA1 = CA1 - N1;
    CartesianPoint Normal1 = N1toCA1.cross(C1 - CA1);
    CartesianPoint N2toCA2 = CA2 - N2;
    CartesianPoint Normal2 = N2toCA2.cross(C2 - CA2);
    
    return cosAngle(Normal1,Normal2);
}

mstreal generalUtilities::cosAngle(CartesianPoint v1, CartesianPoint v2) {
    mstreal dot = v1.dot(v2);
    mstreal v1_mag = v1.norm();
    mstreal v2_mag = v2.norm();
    return dot / (v1_mag * v2_mag);
}

mstreal generalUtilities::avgCosAngleBetweenSegments(const vector<Residue *> &seg1, const vector<Residue *> &seg2) {
    MstUtils::assert(seg1.size() == seg2.size(),"Two segments must have the same number of residues");
    mstreal avg = 0;
    for (int idx = 0; idx < seg1.size(); idx++) {
        avg += cosAngleBetweenNormalVectors(seg1[idx],seg2[idx]);
    }
    return avg / mstreal(seg1.size());
}

bool generalUtilities::contiguousResidues(vector<Residue*> segment, mstreal maxPeptideBond) {
    if (segment.empty()) MstUtils::error("Empty vector of residues was passed","generalUtilities::contiguousResidues");
    for (int i = 0; i < segment.size() - 1; i++) {
        if (!Residue::areBonded(segment[i],segment[i+1],maxPeptideBond)) {
            return false;
        }
    }
    return true;
}

mstreal generalUtilities::bestRMSD(Chain* C1, Chain* C2) {
    if ((C1->residueSize() == 0) || (C2->residueSize() == 0)) MstUtils::error("Both chains must have at least one residue","generalUtilities::bestRMSD");
    
    // Make sure C1 is shorter than C2
    if (C1->residueSize() > C2->residueSize()) {
        Chain* C_temp = C1;
        C1 = C2;
        C2 = C_temp;
    }
    
    // Find the lowest RMSD alignment by sliding C1 over C2
    vector<Residue*> R1 = C1->getResidues();
    vector<Residue*> R2 = C2->getResidues();
    mstreal lowest_rmsd = numeric_limits<double>::max();
    for (int i = 0; i < C2->residueSize() - C1->residueSize() + 1; i++) {
        vector<Residue*> R2_segment = vector<Residue*>(R2.begin()+i,R2.begin()+i+C1->residueSize());
        mstreal new_rmsd = peptideAlignmentNW::getDistance(R1, R2_segment);
        if (new_rmsd < lowest_rmsd) lowest_rmsd = new_rmsd;
    }
    return lowest_rmsd;
}

mstreal volumeCalculator::occupiedVolume(Residue *R) {
    if (protein == nullptr) MstUtils::error("Structure must be set","volumeCalculator::occupiedVolume");
    // find the total volume occupied by atoms within the sphere around the residue alpha carbon
    mstreal occupiedVol = 0.0;
    mstreal maxRadii = vdwRadii::maxRadii();
    Atom* alphaCarbon = R->findAtom("CA");
    if (alphaCarbon == NULL) MstUtils::error("Residue "+R->getChainID()+MstUtils::toString(R->getNum())+" does not have an atom with name 'CA'","volumeCalculator::fracUnoccupiedVolume");
    CartesianPoint CaCoord = alphaCarbon->getCoor();
    
    // Extend the range, since atoms with centers outside of the query sphere could still have some
    // volume within
    vector<int> atomsAroundCenter = ps.getPointsWithin(CaCoord, 0, qR+maxRadii);
    for (int i : atomsAroundCenter) {
        // find the volume of the atom near the center
        Atom* nearAtom = atoms[i];
        CartesianPoint nearAtomCoord = nearAtom->getCoor();
        mstreal r2 = resVDWRadii::getRadii(nearAtom->getResidue()->getName(), nearAtom->getName());
        mstreal distance = alphaCarbon->distance(nearAtom);
        mstreal vdwV = 0.0;
        if (distance >= qR + r2) continue; //atom does not intersect query sphere
        else if (distance < qR - r2) {
            // atom is completely within query sphere
            vdwV = sphereVol(r2);
            occupiedVol += vdwV;
        } else {
            // atom is partially within query sphere
            vdwV = intersectingSphereVol(qR,r2,distance);
            occupiedVol += vdwV;
        }
        
        // find the pairwise overlaps between atom near center and other atoms near the center
        vector<int> atomsAroundNearAtom = ps.getPointsWithin(nearAtomCoord, 0, r2+maxRadii);
        mstreal sharedOccupiedVolume = 0;
        mstreal vdwSharedV = 0.0;
        for (int j : atomsAroundNearAtom) {
            Atom* atomIntersectingNearAtom = atoms[j];
            mstreal r3 = resVDWRadii::getRadii(atomIntersectingNearAtom->getResidue()->getName(), atomIntersectingNearAtom->getName());
            mstreal atomIntersectingNearAtomDistance = nearAtom->distance(atomIntersectingNearAtom);
            if (atomIntersectingNearAtom == nearAtom) continue; // atom intersecting itself
            if (atomIntersectingNearAtomDistance >= r2 + r3) continue; //atom does not intersect query sphere
            mstreal atomIntersectingDistanceToCenter = alphaCarbon->distance(atomIntersectingNearAtom);
            if (atomIntersectingDistanceToCenter >= qR) continue; // intersecting atom mostly outside of sphere
            vdwSharedV = intersectingSphereVol(r2,r3,atomIntersectingNearAtomDistance)/2;
            occupiedVol -= vdwSharedV;
        }
    }
    return occupiedVol;
}

mstreal peptideAlignmentNW::findOptimalAlignment(Structure* _S1, Structure* _S2) {
    // check that these look like normal fused paths
    if ((_S1->chainSize() != 1) || (_S2->chainSize() != 1)) MstUtils::error("One or both of the paths have more than one chain!","peptideAlignmentNW::findOptimalAlignment");
    
    S1 = _S1; //length N ("columns")
    S2 = _S2; //length M ("rows")
    R1 = S1->getResidues();
    R2 = S2->getResidues();
    
    if (verbose) cout << "S1 residue size : " << R1.size() << "\t S2 residue size: " << R2.size() << endl;
    
//    if (S1->residueSize() == S2->residueSize()) {
//        // alignment is trivial
//        best_alignment_distance = getDistance(R1,R2);
//    }
    
    // Resize the tables/alignment vector
    score_table.clear();
    traceback_table.clear();
    score_table.resize(R1.size()+1,vector<mstreal>(R2.size()+1,0.0));
    traceback_table.resize(R1.size(),vector<pair<int,int>>(R2.size(),pair<int,int>(0,0)));
    
    // Fill in the score table (special cases)
    for (int i = 1; i < R1.size()+1; i++) {
        // fill in the leftmost column
        score_table[i][0] = score_table[i-1][0] + gap_penalty;
    }
    for (int j = 1; j < R2.size()+1; j ++) {
        // fill in the topmost row
        score_table[0][j] = score_table[0][j-1] + gap_penalty;
    }

    // Apply the recursive rule to fill the remaining positions
    mstreal gap_chain1, same, gap_chain2;
    for (int i = 1; i < R1.size()+1; i++) {
        for (int j = 1; j < R2.size()+1; j++) {
            // Compute three values corresponding to:
            // 1) gap in chain 1
            gap_chain1 = score_table[i-1][j] + gap_penalty;
            // 2) aligned position
            same = score_table[i-1][j-1] + similarity(R1[i-1],R2[j-1]);
            // 3) gap in chain 2
            gap_chain2 = score_table[i][j-1] + gap_penalty;
            
            // add to score table
            score_table[i][j] = max(gap_chain1,max(same,gap_chain2));
            
            // depending on which value was largest, fill in traceback table
            if ((gap_chain1 > same) && (gap_chain1 > gap_chain2)) {
                traceback_table[i-1][j-1] = pair<int,int>(-1,0);
            } else if ((same > gap_chain1) && (same > gap_chain2)) {
                traceback_table[i-1][j-1] = pair<int,int>(-1,-1);
            } else {
                traceback_table[i-1][j-1] = pair<int,int>(0,-1);
            }
        }
    }
    
    if (verbose) {
        cout.precision(3);
        
        // If verbose, print out score/trackback table
        for (int i = 0; i < R1.size()+1; i++) {
            for (int j = 0; j < R2.size()+1; j++) {
                cout << score_table[i][j] << "\t";
            }
            cout << endl;
        }
        
        for (int i = 0; i < R1.size(); i++) {
            for (int j = 0; j < R2.size(); j++) {
                cout << traceback_table[i][j].first << "," << traceback_table[i][j].second << "\t";
            }
            cout << endl;
        }
    }
    
    // Starting at the (S1_residueSize,S2_residueSize) position in the traceback table, find the path
    // that maximizes the score
    alignment.clear();
    int i = R1.size() - 1; int j = R2.size() - 1;
    while (i >= 0 and j >= 0) {
        pair<int,int>& move = traceback_table[i][j];
        if ((move.first == -1) && (move.second == -1)) {
            if (verbose) cout << "\ti: " << i << "\tj: " << j << endl;
            alignment.push_front(pair<int,int>(i,j));
        }
        i += move.first;
        j += move.second;
    }
    
    // Get residues using alignment
    R1_aligned.clear(); R2_aligned.clear();
    for (pair<int,int>& alignment_pos : alignment) {
        if (verbose) cout << alignment_pos.first << " " << alignment_pos.second << endl;
        R1_aligned.push_back(R1[alignment_pos.first]);
        R2_aligned.push_back(R2[alignment_pos.second]);
    }
    
    // Get RMSD
    best_alignment_distance = getDistance(R1_aligned,R2_aligned);
    
    if (verbose) cout << "distance: " << best_alignment_distance << endl;
    
    return best_alignment_distance;
}

mstreal peptideAlignmentNW::getDistance(vector<Residue*> R1, vector<Residue*> R2) {
    vector<Atom*> A1; vector<Atom*> A2;
    for (int i = 0; i < R1.size(); i++) {
        Residue* R1_i = R1[i];
        Residue* R2_i = R2[i];
        if (!RotamerLibrary::hasFullBackbone(R1_i) || !RotamerLibrary::hasFullBackbone(R2_i)) MstUtils::error("Residue from peptide is missing backbone atoms","peptideAlignmentNW::getRMSD");
        vector<Atom*> R1_i_bb_atoms = RotamerLibrary::getBackbone(R1_i);
        vector<Atom*> R2_i_bb_atoms = RotamerLibrary::getBackbone(R2_i);
        A1.insert(A1.end(),R1_i_bb_atoms.begin(),R1_i_bb_atoms.end());
        A2.insert(A2.end(),R2_i_bb_atoms.begin(),R2_i_bb_atoms.end());
    }
    RMSDCalculator rmsd_calc;
    return rmsd_calc.rmsd(A1,A2);
}

mstreal peptideAlignmentNW::similarity(Residue* R1, Residue* R2) {
    if (!RotamerLibrary::hasFullBackbone(R1) || !RotamerLibrary::hasFullBackbone(R2)) MstUtils::error("Residue from peptide is missing backbone atoms","peptideAlignmentNW::similarity");
    vector<Atom*> R1_bb_atoms = RotamerLibrary::getBackbone(R1);
    vector<Atom*> R2_bb_atoms = RotamerLibrary::getBackbone(R2);
    RMSDCalculator rmsd_calc;
    return -rmsd_calc.rmsd(R1_bb_atoms,R2_bb_atoms);
}

/* ----------- Miscellaneous useful functions -------------- */

int getTargetResidueIndex(string seedName) {
    regex re("^[A-z0-9]{4}_[A-z](\\d+)_");
    smatch match;
    string fileName = MstSystemExtension::fileName(seedName);
    if (std::regex_search(fileName, match, re) && match.size() > 1) {
        return atoi(match.str(1).c_str());
    }
    return -1;
}

pair<string, int> getTargetResidueCode(string seedName) {
    regex re("([A-z])(\\d+)");
    smatch match;
    string fileName = MstSystemExtension::fileName(seedName);
    string central_residue = MstUtils::split(fileName,"-")[1];
    if (std::regex_search(central_residue, match, re) && match.size() > 1) {
        return make_pair(match.str(1), atoi(match.str(2).c_str()));
    }
    return make_pair("", -1);
}

vector<string> splitString(string s, const string &delim) {
    size_t pos = 0;
    string token;
    vector<string> result;
    while ((pos = s.find(delim)) != string::npos) {
        token = s.substr(0, pos);
        result.push_back(token);
        s.erase(0, pos + delim.length());
    }
    result.push_back(s);
    return result;
}


// ========= Utilities from structgen

// tests whether atom is backbone
// no integer second argument in this function
bool isBackboneHeavy(Atom& atom) {
    return (atom.getName() == "N" || atom.getName() == "NT" || atom.getName() == "CA" || atom.getName() == "C" || atom.getName() == "CT" || atom.getName() == "O");
}

// determines if a residue contains all of its backbone atoms
bool hasBackbone(Residue& res, bool requireOxygen) {
    return ((res.atomExists("N") || res.atomExists("NT")) &&
            res.atomExists("CA") &&
            (res.atomExists("C") || res.atomExists("CT")) &&
            (res.atomExists("O") || !requireOxygen));
}

// distance between the C termini of res1 and the N termini of res2
double bondDistance(Residue& res1, Residue& res2) {
    Atom* c = res1.findAtom("C", false);
    if (c == NULL) {
        c = res1.findAtom("CT", true);
    }
    Atom* n = res2.findAtom("N", false);
    if (n == NULL) {
        n = res2.findAtom("NT", true);
    }
    return c->distance(n);
}

// is res1 close enough to be bonded to res2?
// order matters here
bool isBonded(Residue& res1, Residue& res2, double maxPepBond) {
    return (bondDistance(res1, res2) <= maxPepBond);
}

// removes the side chains from a structure
void removeSideChains(Structure& structure) {
    vector<Residue*> residues = structure.getResidues();
    for (auto res = residues.begin(); res != residues.end(); res++) {
        int j = 0;
        while (j < (*res)->atomSize()) {
            if (!isBackboneHeavy((*res)->getAtom(j))) {
                (*res)->deleteAtom(j);
            } else {
                j++;
            }
        }
    }
}

// returns a new AtomPointerVector without the side chains
AtomPointerVector backboneAtoms(AtomPointerVector& apv) {
    AtomPointerVector mainChain;
    for(int i = 0; i < apv.size(); i++) {
        if (isBackboneHeavy(*apv[i])) {
            mainChain.push_back(apv[i]);
        }
    }
    return mainChain;
}

// removes the hydrogens from a structure
void removeHydrogens(Structure& s) {
    vector<Residue*> residues = s.getResidues();
    for (auto res = residues.begin(); res != residues.end(); res++) {
        int j = 0;
        while (j < (*res)->atomSize()) {
            string atomName = (*res)->getAtom(j).getName();
            if (atomName[0] == 'H' || (atomName.length() > 1 && atomName[1] == 'H')) {
                (*res)->deleteAtom(j);
            } else {
                j++;
            }
        }
    }
}

// an alphabet for naming chains
string generateChainID(int i) {
    string alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    stringstream ss;
    ss << alphabet[i % alphabet.size()];
    return ss.str();
}

// combines a vector of structures into a single structure where each former structure is now a chain
// will always renumber if split is false
Structure combineStructures(const vector<Structure>& structs, bool split, bool renumber) {
    Structure comb;
    for (int i = 0; i < structs.size(); i++) {
        Structure s = structs[i];
        if (split) { // splits by existing chains and not by connectivity
            for (int j = 0; j < s.chainSize(); j++) {
                Chain* c = comb.appendChain(generateChainID(comb.chainSize()));
                c->setSegID(c->getID());
                c->appendResidueCopies(s[j].getResidues());
            }
        } else {
            Chain* c = comb.appendChain(generateChainID(comb.chainSize()));
            c->setSegID(c->getID());
            vector<Residue*> residues = s.getResidues();
            for (int j = 0; j < residues.size(); j++) {
                c->appendResidue(new Residue(*residues[j]));
                Residue* res = &c->getResidue(j);
                res->setNum(j + 1);
            }
        }
    }
    if (renumber) {
        comb.renumber();
    }
    return comb;
}

// splits a structure into a vector of structures based on the chains
vector<Structure> splitByChains(const Structure& s) {
    vector<Structure> structs;
    for (int i = 0; i < s.chainSize(); i++) {
        Chain c = s[i];
        c.setID("A"); // is this nesecary
        Structure si(c);
        si.renumber();
        structs.push_back(si);
    }
    return structs;
}


string residueID(Residue& res, string sep, bool full) {
    string segID = "";
    if (full && res.getParent()->getSegID() != "") {
        segID += res.getParent()->getSegID() + sep;
    }
    
    string icode = "";
    if (full && res.getIcode() != ' ') {
        icode = sep + res.getIcode();
    }
    return segID + res.getChainID() + sep + MstUtils::toString(res.getNum()) + icode;
}


// writes residues from a vector out
void writeResidues(ostream& os, const vector<Residue*>& residues) {
    for (int i = 0; i < residues.size(); i++) {
        os << *residues[i] << endl;
    }
}

// writes residues from a vector out
void writeResidues(ostream& os, const vector<vector<Residue*> >& residues) {
    for (int i = 0; i < residues.size(); i++) {
        for (int j = 0; j < residues[i].size(); j++) {
            os << *residues[i][j];
            if (j < residues[i].size() - 1) {
                cout << " ";
            }
        }
        os << endl;
    }
}

// sorts a vector of residues by their index
vector<Residue*> sortByResidueIndex(vector<Residue*>& residues) {
    map<int, Residue*> indices;
    for (int i = 0; i < residues.size(); i++) {
        indices[residues[i]->getResidueIndex()] = residues[i];
    }
    vector<Residue*> sortedResidues;
    for (auto it = indices.begin(); it != indices.end(); it++) {
        sortedResidues.push_back(it->second);
    }
    return sortedResidues;
}


vector<Residue*> sortByResidueIndex(set<Residue*>& residues) {
    vector<Residue*> tmp(residues.begin(), residues.end());
    return sortByResidueIndex(tmp);
}


// map of Residue -> index within structure
// fromResidue means to get the index from the residue rather than from its position within the vector of residues
map<Residue*, int> residueIndexMapping(vector<Residue*>& residues, bool fromResidue) {
    map<Residue*, int> resMapping;
    for (int i = 0; i < residues.size(); i++) {
        if (fromResidue) {
            resMapping[residues[i]] = residues[i]->getResidueIndex();
        } else {
            resMapping[residues[i]] = i;
        }
    }
    return resMapping;
}

// converts a vector of Residues into an AtomPointerVector
AtomPointerVector residuesToAtoms(vector<Residue*>& residues) {
    AtomPointerVector apv;
    for (int i = 0; i < residues.size(); i++) {
        AtomPointerVector resAtoms = residues[i]->getAtoms();
        apv.insert(apv.end(), resAtoms.begin(), resAtoms.end());
    }
    return apv;
}

// calculates the sequence frequencies from a set of sequences
Matrix seqFreqs(const vector<Sequence>& seqs, double pseudocount, bool normalize, bool nonStd) {
    MstUtils::assert(seqs.size() > 0, "Error in seqFreqs(vector<Sequence>& seqs, bool normalize, double pseudocounts). An empty vector of sequences was passed.");
    int numRows = ((nonStd) ? 21 : 20);
    Matrix psuedocounts(numRows, seqs[0].length());
    for (int i = 0; i < 20; i++) { // don't add pseudocount for 21st row if there is one
        for (int j = 0; j < seqs[0].length(); j++) {
            psuedocounts(i, j) = pseudocount;
        }
    }
    
    return seqFreqs(seqs, psuedocounts, normalize);
}


Matrix seqFreqs(const vector<Sequence>& seqs, Matrix& pseudocounts, bool normalize) {
    MstUtils::assert(seqs.size() > 0, "Error in seqFreqs(vector<Sequence>& seqs, bool normalize, Matrix& pseudocounts). An empty vector of sequences was passed.");
    int seqLen = seqs[0].length();
    MstUtils::assert(pseudocounts.numCols() == seqLen && (pseudocounts.numRows() == 20 || pseudocounts.numRows() == 21), "Error in seqFreqs(vector<Sequence>& seqs, bool normalize, Matrix& pseudocounts). Matrix of pseudocounts must have 20 rows and the number of columns must be the same as the sequence length.");
    
    Matrix freqs(pseudocounts); // 20 x L matrix where L is the length of the sequences
    
    for (int i = 0; i < seqs.size(); i++) {
        Sequence seq = seqs[i];
        MstUtils::assert(seq.length() == seqLen, "Sequences not the same length."); // better place to do this? perhaps alignment class
        
        for (int j = 0; j < seqLen; j++) {
            string a = seq.getResidue(j);
            int aaIdx = SeqTools::aaToIdx(a);
            aaIdx = (aaIdx >= 20 && freqs.numRows() > 20) ? 20 : aaIdx; // if aa index is over 19, it is non-standard, put it in the 20th row if there is a 20th row
            if (aaIdx < freqs.numRows()) { // make sure it is a standard amino acid
                freqs(aaIdx, j) += 1.0;
            }
        }
    }
    
    if (normalize) { // normalize the frequency matrix (column-wise)
        for (int c = 0; c < seqLen; c++) {
            // find sums of columns
            double colSum = 0.0;
            for (int r = 0; r < freqs.numRows(); r++) {
                colSum += freqs(r, c);
            }
            if (colSum > 0) { // divide by sum of column
                for (int r = 0; r < freqs.numRows(); r++) {
                    freqs(r, c) /= colSum;
                }
            }
        }
    }
    
    return freqs;
}


// remove chains below a certain size (e.g. 2 residues)
// option for dealing non standard AAs
// incorporate TargetStruct redundancy here
void cleanStructure(Structure& s, Structure& cleaned, bool reassignChains, bool renumber, int minSegLen) {
    cleaned.reset();
    vector<Residue*> residues = s.getResidues();
    
    vector<Residue*> passingResidues;
    // check that all backbone atoms are present
    for (int i = 0; i < residues.size(); i++) {
        Residue* res = residues[i];
        if (res->atomSize() >= 4 && hasBackbone(*res)) {
            //previously residues were not added if they were named 'UNK' or '-'.
//            res_t idx = SeqTools::aaToIdx(res->getName());
//            if (idx != SeqTools::unknownIdx() && idx != SeqTools::gapIdx()) {
//                passingResidues.push_back(res);
//            }
            passingResidues.push_back(res);
        }
    }
    cleaned = Structure(passingResidues);
    if (minSegLen > 1) {
        vector<Residue*> cleanedResidues = cleaned.getResidues();
        passingResidues.clear();
        vector<Residue*> segResidues;
        for (int i = 0; i < cleanedResidues.size(); i++) {
            segResidues.push_back(cleanedResidues[i]);
            if ((i == cleanedResidues.size() - 1) || !isBonded(*cleanedResidues[i], *cleanedResidues[i + 1])) {
                if (segResidues.size() >= minSegLen) {
                    passingResidues.insert(passingResidues.end(), segResidues.begin(), segResidues.end());
                }
                segResidues.clear();
            }
        }
        cleaned = Structure(passingResidues);
    }
    if (reassignChains) { // reassign chains based on peptide bond distances
        cleaned = cleaned.reassignChainsByConnectivity();
    }
    
    // renumber residues within chains
    if (renumber) {
        cleaned.renumber();
    }
}

Structure cleanStructure(Structure& s, bool reassignChains, bool renumber, int minSegLen) {
    Structure cleaned;
    cleanStructure(s, cleaned, reassignChains, renumber, minSegLen);
    return cleaned;
}


// vector of chain lengths for a structure
vector<int> chainLengths(Structure& s) {
    vector<int> v;
    for (int i = 0; i < s.chainSize(); i++) {
        v.push_back(s[i].residueSize());
    }
    return v;
}

// calulate density
double calculateDensity(Structure& s) {
    AtomPointerVector apv = s.getAtoms();
    double rad = apv.boundingSphereRadiusCent();
    return s.residueSize() / ((4.0 / 3.0) * M_PI * rad * rad * rad);
}

// calculate pairwise RMSD
void pairwiseRMSD(vector<AtomPointerVector>& apvs, unordered_map<int, unordered_map<int, double>>& rmsds, double rmsdCut) {
    int i, j;
    double rmsd;
    RMSDCalculator rmsdCalc;
    bool suc;
    for (i = 0; i < apvs.size() - 1; i++) {
        for (j = i + 1; j < apvs.size(); j++) {
            rmsdCalc.bestRMSD(apvs[i], apvs[j], false, &suc);
            rmsd = rmsdCalc.lastRMSD();
            if (suc && rmsd < rmsdCut) {
                rmsds[i][j] = rmsd;
                rmsds[j][i] = rmsd;
            }
        }
    }
}


set<pair<Atom*, Atom*>> findClashes(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& qAPV, double ratio) {
    set<pair<Atom*, Atom*>> clashes;
    
    // this should be passed in so that user can specify these - e.g. could include side chain atoms as well
    double dist = vdwRadii::maxSumRadii() * ratio;
    for (int i = 0; i < qAPV.size(); i++) {
        Atom* ai = qAPV[i];
        Residue* resi = ai->getParent();
        vector<int> tags = ps.getPointsWithin(ai->getCoor(), 0, dist, true);
        for (int j = 0; j < tags.size(); j++) {
            Atom* aj = psAPV[tags[j]];
            Residue* resj = aj->getParent();
            if (resi != resj && // not in the same residue
                (resi->previousResidue() != resj || !isBonded(*resj, *resi)) && // not in the prev bonded residue
                (resi->nextResidue() != resj || !isBonded(*resi, *resj)) ) { // not in next bonded residue
                if (vdwRadii::clash(*ai, *aj, ratio)) {
                    clashes.insert(orderedPair(ai, aj));
                }
            }
        }
    }
    return clashes;
}


// simple clash function
// does not take into account atom sizes, ect
// no ability to exclude certain atoms or residues
bool isClash(ProximitySearch& ps, AtomPointerVector& queryAPV, double clashDist, int maxNumClashes) {
    int numClashes = 0;
    for (int i = 0; i < queryAPV.size(); i++) {
        vector<int> tags = ps.getPointsWithin(queryAPV[i]->getCoor(), 0, clashDist);
        numClashes += tags.size();
        if (numClashes > maxNumClashes) {
            return true;
        }
    }
    return false;
}

// determines whether there a clash betwen an AtomPointerVector (apv) and another structure represented by a ProximitySearch object (ps)
// ps is a ProximitySearchObject
// apv is an AtomPointerVector
// psExclude is the set of atom indices in psAPV (and ps) to exlude from clashes
// ratio of sum of vdw radii (default = 0.7)
// maxNumClasshes is the maximum number of allowed clashes between queryAPV and psAPV (default = 0)
// raddii comes from /home/ironfs/scratch/grigoryanlab/cmack2357/minCover/lists
bool isClash(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& queryAPV, set<pair<Residue*, Residue*> >& exclude, double ratio, int maxNumClashes) {
    int numClashes = numClash(ps, psAPV, queryAPV, exclude, ratio, maxNumClashes);
    if (numClashes > maxNumClashes) return true;
    else return false;
}

// calls the above function without and excluded residues
bool isClash(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& queryAPV, double ratio, int maxNumClashes) {
    set<pair<Residue*, Residue*> > exclude;
    return isClash(ps, psAPV, queryAPV, exclude, ratio, maxNumClashes);
}

int numClash(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& queryAPV, set<pair<Residue*, Residue*> >& exclude, double ratio, int maxNumClashes) {
    int numClashes = 0;
    
    // this should be passed in so that user can specify these - e.g. could include side chain atoms as well
    double maxDist = vdwRadii::maxSumRadii() * ratio;
    
    for (int qi = 0; qi < queryAPV.size(); qi++) {
        Residue* qRes = queryAPV[qi]->getParent();
        vector<int> tags = ps.getPointsWithin(queryAPV[qi]->getCoor(), 0, maxDist); //, true);
        for (int j = 0; j < tags.size(); j++) {
            int pi = tags[j];
            Residue* psRes = psAPV[pi]->getParent();
            //check if either direction of residue clash is excluded
            pair<Residue*,Residue*> psRes_qRes(psRes, qRes);
            pair<Residue*,Residue*> qRes_psRes(qRes, psRes);
            if ((exclude.count(psRes_qRes) < 1) && (exclude.count(qRes_psRes) < 1)) {
                if (abs(qRes->getResidueIndex() - psRes->getResidueIndex()) <= 1) continue;
                //if both directions were not in excluded, check if counts as clash
                if (vdwRadii::clash(*psAPV[pi], *queryAPV[qi], ratio)) {
                    numClashes ++;
                    //add to excluded so it is not counted again later
                    exclude.insert(make_pair(psRes, qRes));
                    exclude.insert(make_pair(qRes, psRes));
                }
                if ((maxNumClashes >= 0) & (numClashes > maxNumClashes)) {
                    return numClashes;
                }
            }
        }
    }
    return numClashes;
}

//// Added by venkats 1/28/20
//// Same as above, but operates on a single structure
//bool isClashSingleStructure(ProximitySearch& ps, AtomPointerVector& psAPV, double ratio, int maxNumClashes) {
//    int numClashes = numClashSingleStructure(ps, psAPV, ratio, maxNumClashes);
//    if (numClashes >= maxNumClashes) return true;
//    else return false;
//}
//
//int numClashSingleStructure(ProximitySearch& ps, AtomPointerVector& psAPV, double ratio, int maxNumClashes) {
//    int numClashes = 0;
//
//    // this should be passed in so that user can specify these - e.g. could include side chain atoms as well
//    double maxDist = vdwRadii::maxSumRadii() * ratio;
//
//    for (int qi = 0; qi < psAPV.size(); qi++) {
//        Residue* qRes = psAPV[qi]->getParent();
//        vector<int> tags = ps.getPointsWithin(psAPV[qi]->getCoor(), 0, maxDist); //, true);
//        for (int j = 0; j < tags.size(); j++) {
//            int pi = tags[j];
//            Residue* psRes = psAPV[pi]->getParent();
//            if (abs(qRes->getResidueIndex() - psRes->getResidueIndex()) <= 1)
//                continue;
//            if (vdwRadii::clash(*psAPV[pi], *psAPV[qi], ratio)) {
//                numClashes ++;
//            }
//            if ((maxNumClashes >= 0) & (numClashes > maxNumClashes)) {
//                return numClashes;
//            }
//        }
//    }
//    return numClashes;
//}

vector<int> effectiveSegLengths(Structure& s, vector<Residue*> subset, int minGapCont, int minSegLen) {
    vector<Residue*> residues = s.getResidues();
    vector<int> segLens;
    int segLen = 1;
    for (int i = 1; i < subset.size(); i++) {
        Residue* prevRes = residues[i - 1];
        Residue* curRes = residues[i];
        int prevResi = s.getResidueIndex(prevRes);
        int curResi = s.getResidueIndex(curRes);
        
        int gap = prevResi - (curResi + 1);
        bool connected = true;
        int numBonds = 0;
        if (gap == 0 || gap < minGapCont) { // the difference is less than the minimum difference to called separate - now we just need to see if they are connected
            for (int j = prevResi; j < curResi; j++) {
                Residue* res1 = residues[j];
                Residue* res2 = residues[j + 1];
                if (!isBonded(*res1, *res2)) {
                    connected = false;
                    numBonds = 0;
                    break;
                }
                numBonds++;
            }
        }
        segLen += numBonds;
        if (!connected || (i == subset.size() - 1)) {
            segLens.push_back(max(segLen, minSegLen));
            segLen = 1;
        }
    }
    return segLens;
}

Structure residuesToStructure(vector<Residue*>& residues, double maxPeptideBond, int startChainIdx) {
    
    Structure s;
    if (residues.size() == 0) {
        return s;
    }
    
    int chainIdx = startChainIdx;
    Chain* chain = s.appendChain(generateChainID(chainIdx));
    
    for (int i  = 0; i < residues.size() - 1; i++) {
        Atom* atomC = residues[i]->findAtom("C", true);
        Atom* atomN = residues[i + 1]->findAtom("N", true);
        double dist = atomC->distance(atomN);
        chain->appendResidue(new Residue(*residues[i]));
        if (dist > maxPeptideBond) {
            chainIdx++;
            chain = s.appendChain(generateChainID(chainIdx));
        }
    }
    chain->appendResidue(new Residue(*residues[residues.size() - 1]));
    s.renumber();
    return s;
}

// writes the contents of a contactList
void writeContactList(ostream& os, contactList& cl) {
    for (int i = 0; i < cl.size(); i++) {
        os << *cl.residueA(i) << " " << *cl.residueB(i) << " " << cl.degree(i) << endl;
    }
}


// get the backbone backbone contacts for a structure
// queryResidues should be a subset of the residues in s
contactList vdwBackboneContacts(Structure& s, double lb, double ub, vector<Residue*> queryResidues) {
    
    AtomPointerVector apv = s.getAtoms();
    AtomPointerVector mainChain = backboneAtoms(apv);
    ProximitySearch ps(mainChain, 40.0);
    
    vector<Residue*> residues;
    if (queryResidues.size() == 0) {
        residues = s.getResidues();
    } else {
        residues = queryResidues;
        // need to check the residues are in s
    }
    
    map<pair<Residue*, Residue*>, double> seen;
    for (int i = 0; i < residues.size(); i++) {
        Residue* qRes = residues[i];
        AtomPointerVector qAtoms = qRes->getAtoms();
        AtomPointerVector bbQAtoms = backboneAtoms(qAtoms);
        
        for (int j = 0; j < bbQAtoms.size(); j++) {
            Atom* qAtom = bbQAtoms[j];
            vector<int> points = ps.getPointsWithin(qAtom, 1.0, 4.0); // I guess it should also be above the clash filter dist
            Residue* refRes = qAtom->getParent();
            
            for (int k = 0; k < points.size(); k++) {
                Atom* refAtom = mainChain[points[k]];
                Residue* refRes = refAtom->getParent();
                if (refRes != qRes) {
                    int refResi = s.getResidueIndex(refRes);
                    if ((refResi != i - 1 || !isBonded(*refRes, *qRes)) && (refResi != i + 1 || !isBonded(*qRes, *refRes))) {
                        pair<Residue*, Residue*> p = orderedPair(qRes, refRes);
                        double d = qAtom->distance(refAtom);
                        bool contacting = resVDWRadii::contact(*qAtom, *refAtom, lb, ub);
                        if (contacting && (seen.count(p) == 0 || d < seen[p])) {
                            seen[p] = d;
                        }
                    }
                }
            }
        }
    }
    
    contactList cl;
    for (auto it = seen.begin(); it != seen.end(); it++) {
        cl.addContact(it->first.first, it->first.second, it->second);
    }
    
    return cl;
}

// get the backbone backbone contacts for a structure
// queryResidues should be a subset of the residues in s
contactList backboneContacts(Structure& s, double dist, vector<Residue*> queryResidues) {
    
    AtomPointerVector apv = s.getAtoms();
    AtomPointerVector mainChain = backboneAtoms(apv);
    ProximitySearch ps(mainChain, 40.0);
    
    vector<Residue*> residues;
    if (queryResidues.size() == 0) {
        residues = s.getResidues();
    } else {
        residues = queryResidues;
        // need to check the residues are in s
    }
    
    map<pair<Residue*, Residue*>, double> seen;
    for (int i = 0; i < residues.size(); i++) {
        Residue* qRes = residues[i];
        AtomPointerVector qAtoms = qRes->getAtoms();
        AtomPointerVector bbQAtoms = backboneAtoms(qAtoms);
        
        for (int j = 0; j < bbQAtoms.size(); j++) {
            Atom* qAtom = bbQAtoms[j];
            vector<int> points = ps.getPointsWithin(qAtom, 1.0, dist); // I guess it should also be above the clash filter dist
            Residue* refRes = qAtom->getParent();
            
            for (int k = 0; k < points.size(); k++) {
                Atom* refAtom = mainChain[points[k]];
                Residue* refRes = refAtom->getParent();
                if (refRes != qRes) {
                    int refResi = s.getResidueIndex(refRes);
                    if ((refResi != i - 1 || !isBonded(*refRes, *qRes)) && (refResi != i + 1 || !isBonded(*qRes, *refRes))) {
                        pair<Residue*, Residue*> p = orderedPair(qRes, refRes);
                        double d = qAtom->distance(refAtom);
                        if (seen.count(p) == 0 || d < seen[p]) {
                            seen[p] = d;
                        }
                    }
                }
            }
        }
    }
    
    contactList cl;
    for (auto it = seen.begin(); it != seen.end(); it++) {
        cl.addContact(it->first.first, it->first.second, it->second);
    }
    
    return cl;
}

// union of contact list
// does not allow duplicates
// orders the contacts - perhaps give an option to not do this
contactList contactListUnion(const vector<contactList>& contLists, bool order) {
    set<pair<Residue*, Residue*> > seen;
    contactList combined;
    for (int i = 0; i < contLists.size(); i++) {
        contactList cl = contLists.at(i);
        for (int j = 0; j < cl.size(); j++) {
            pair<Residue*, Residue*> cont;
            if (order) cont = orderedPair(cl.residueA(j), cl.residueB(j));
            else cont = make_pair(cl.residueA(j), cl.residueB(j));
            if (seen.count(cont) == 0) {
                combined.addContact(cont.first, cont.second, cl.degree(j));
                seen.insert(cont);
            }
        }
    }
    return combined;
}

