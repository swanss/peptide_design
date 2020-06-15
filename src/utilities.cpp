//
//  utilities.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 6/14/19.
//

#include "utilities.h"
#include "mstsystem_exts.h"
#include "vdwRadii.h"
#include "mstlinalg.h"
#include <regex>

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
    regex re("^[A-z0-9_]-([A-z])(\\d+)-");
    smatch match;
    string fileName = MstSystemExtension::fileName(seedName);
    if (std::regex_search(fileName, match, re) && match.size() > 1) {
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
        if (res->atomSize() >= 4 && hasBackbone(*res) && res->getName().length() == 3) {
            res_t idx = SeqTools::aaToIdx(res->getName());
            if (idx != SeqTools::unknownIdx() && idx != SeqTools::gapIdx()) {
                passingResidues.push_back(res);
            }
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
    int numClashes = 0;

    // this should be passed in so that user can specify these - e.g. could include side chain atoms as well
    double maxDist = vdwRadii::maxSumRadii() * ratio;

    for (int qi = 0; qi < queryAPV.size(); qi++) {
        Residue* qRes = queryAPV[qi]->getParent();
        vector<int> tags = ps.getPointsWithin(queryAPV[qi]->getCoor(), 0, maxDist); //, true);
        for (int j = 0; j < tags.size(); j++) {
            int pi = tags[j];
            Residue* psRes = psAPV[pi]->getParent();
            if (exclude.count(make_pair(psRes, qRes)) == 0) {
                if (vdwRadii::clash(*psAPV[pi], *queryAPV[qi], ratio)) {
                    numClashes ++;
                }
                if (numClashes > maxNumClashes) {
                    return true;
                }
            }
        }
    }
    return false;
}

// Added by venkats 1/28/20
// Same as above, but operates on a single structure
bool isClashSingleStructure(ProximitySearch& ps, AtomPointerVector& psAPV, double ratio, int maxNumClashes) {
    int numClashes = 0;

    // this should be passed in so that user can specify these - e.g. could include side chain atoms as well
    double maxDist = vdwRadii::maxSumRadii() * ratio;

    for (int qi = 0; qi < psAPV.size(); qi++) {
        Residue* qRes = psAPV[qi]->getParent();
        vector<int> tags = ps.getPointsWithin(psAPV[qi]->getCoor(), 0, maxDist); //, true);
        for (int j = 0; j < tags.size(); j++) {
            int pi = tags[j];
            Residue* psRes = psAPV[pi]->getParent();
            if (abs(qRes->getResidueIndex() - psRes->getResidueIndex()) <= 1)
                continue;
            if (vdwRadii::clash(*psAPV[pi], *psAPV[qi], ratio)) {
                numClashes ++;
            }
            if (numClashes > maxNumClashes) {
                return true;
            }
        }
    }
    return false;
}
// calls the above function without and excluded residues
bool isClash(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& queryAPV, double ratio, int maxNumClashes) {
    set<pair<Residue*, Residue*> > exclude;
    return isClash(ps, psAPV, queryAPV, exclude, ratio, maxNumClashes);
}

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
contactList contactListUnion(const vector<contactList>& contLists) {
    set<pair<Residue*, Residue*> > seen;
    contactList combined;
    for (int i = 0; i < contLists.size(); i++) {
        contactList cl = contLists.at(i);
        for (int j = 0; j < cl.size(); j++) {
            pair<Residue*, Residue*> cont = orderedPair(cl.residueA(j), cl.residueB(j));
            if (seen.count(cont) == 0) {
                combined.addContact(cont.first, cont.second, cl.degree(j));
                seen.insert(cont);
            }
        }
    }
    return combined;
}

