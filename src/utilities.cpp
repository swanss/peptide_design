//
//  utilities.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 6/14/19.
//

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

///* --------- StructuresBinaryFile --------- */
//Structure *StructuresBinaryFile::next() {
//  MstUtils::assert(readMode, "next not supported in write mode");
//  Structure* S = new Structure();
//  S->readData(fs);
//  return S;
//}
//
//Structure * StructuresBinaryFile::getStructureNamed(string name) {
//  MstUtils::assert(readMode, "getStructureNamed not supported in write mode");
//  if (_filePositions.size() == 0) {
//    scanFilePositions();
//  }
//  if (_filePositions.count(name) == 0) {
//    cout << "Structure doesn't exist" << endl;
//    return NULL;
//  }
//  fs.seekg(_filePositions[name], fs.beg);
//  Structure *S = new Structure();
//  S->readData(fs);
//  return S;
//}
//
//bool StructuresBinaryFile::hasNext() {
//  MstUtils::assert(readMode, "hasNext not supported in write mode");
//  return fs.peek() != EOF;
//}
//
//void StructuresBinaryFile::skip() {
//  MstUtils::assert(readMode, "skip not supported in write mode");
//  Structure *s = new Structure();
//  s->readData(fs);
//  delete s;
//}
//
//void StructuresBinaryFile::jumpToStructureIndex(int idx) {
//  MstUtils::assert(readMode, "jumpToStructureIndex not supported in write mode");
//  if (idx < 0)
//    fs.seekg(0, fs.beg);
//  else if (idx >= _structureNames.size())
//    fs.seekg(0, fs.end);
//  else {
//    string name = _structureNames[idx];
//    fs.clear();
//    fs.seekg(_filePositions[name], fs.beg);
//  }
//}
//
//void StructuresBinaryFile::reset() {
//  MstUtils::assert(readMode, "reset not supported in write mode");
//  cout << "Resetting" << endl;
//  fs.clear(); // this is necessary in case ifstream doesn't clear eofbit
//  fs.seekg(0, fs.beg);
//}
//
//void StructuresBinaryFile::appendStructure(Structure *s) {
//  MstUtils::assert(!readMode, "appendStructure not supported in read mode");
//  s->writeData(fs);
//}
//
//void StructuresBinaryFile::scanFilePositions() {
//  //MstUtils::assert(readMode, "scanFilePositions not supported in write mode");
//  if (!_structureNames.empty())
//    return;
//  cout << "Scanning file positions..." << endl;
//  fs.seekg(0, fs.beg);
//  _structureNames.clear();
//  while (fs.peek() != EOF) {
//    Structure *S = new Structure();
//    long pos = fs.tellg();
//    S->readData(fs);
//    if (pos < 0) {
//      cout << pos << endl;
//    }
//    _filePositions[S->getName()] = pos;
//    _structureNames.push_back(S->getName());
//    delete S;
//  }
//  reset();
//  cout << "Done scanning file" << endl;
//}
//
//void StructuresBinaryFile::openFileStream(string filePath) {
//  if (readMode)
//    MstUtils::openFile(fs, filePath, fstream::in | fstream::binary, "StructuresBinaryFile::openFileStream");
//  else
//    MstUtils::openFile(fs, filePath, fstream::out | fstream::binary, "StructuresBinaryFile::openFileStream");
//}

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
