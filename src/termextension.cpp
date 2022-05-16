#include "termextension.h"

using namespace MST;

/* --------- seedTERM --------- */
seedTERM::seedTERM() {
    // compiler error thrown unless I include this
}

seedTERM::seedTERM(TermExtension* FragmenterObj, vector<Residue*> _interactingRes, int _flankRes, bool _search, mstreal max_rmsd, bool discard_match_data) : max_RMSD_cutoff(max_rmsd), flankRes(_flankRes), search(_search), interactingRes(_interactingRes) {
    MstTimer timer;
    
    // set up the fragment
    parent = FragmenterObj;
    cenRes = interactingRes[0];
    cenResName = cenRes->getChainID() + MstUtils::toString(cenRes->getNum());
    num_extended_fragments = 0;
    num_seed_residues = 0;
    
    //generate the structure
    if (parent->params.verbose) cout << "Creating TERM centered around residues: " << MstUtils::vecPtrToString(interactingRes," ") << endl;
    cenResIdx = TERMUtils::selectTERM(interactingRes, fragmentStructure, flankRes, &(allResIdx))[0];
    
    //get the complexity
    effective_degrees_of_freedom = generalUtilities::fragDegreesOfFreedom(allResIdx, *parent->getTarget(), 15);
    rmsd_adjust_factor = 1.0/sqrt(fragmentStructure.atomSize()/effective_degrees_of_freedom);
    
    //generate a name that is unique to the fragment
    //NOTE: this name uses the residue numbers from the parent structure
    setName();
    
    // find the RMSD, regardless of whether fragment will be searched
    adaptive_RMSD_cutoff = RMSDCalculator::rmsdCutoff(allResIdx, *parent->getTarget(), max_RMSD_cutoff, 15); // NOTE: the residues in the same chain of the original structure are not by default assumed to be independent
    RMSD_cutoff = FragmenterObj->getAdaptiveRMSD() ? adaptive_RMSD_cutoff : max_RMSD_cutoff;
    
    if (search) searchForMatches(discard_match_data);
    else construction_time = -1.0;
    //store the secondary structure in the resname field of the fragmentStructure residues;
    vector<Residue*> fragmentResidues = fragmentStructure.getResidues();
    for (int i = 0; i < allResIdx.size(); i++) {
        Residue* R_original = &(parent->target_structure->getResidue(allResIdx[i]));
        Residue* R_fragmentCopy = fragmentResidues[i];
        R_fragmentCopy->setName(parent->sec_structure.classifyResidue(R_original));
    }
};

seedTERM::~seedTERM() {
    for (int i = 0; i < seeds.size(); i++) {
        delete seeds[i].extended_fragment;
    }
}

void seedTERM::setName() {
    stringstream ss;
    // parent structure name
    ss << parent->getName() << "-";
    
    // all contacting residue chain AND number
    for (int i = 0; i < interactingRes.size(); i++) {
        ss << interactingRes[i]->getChainID() << MstUtils::toString(interactingRes[i]->getNum());
        if (i == 0) ss << "-"; //central residue is always first, separate from the rest with the hyphen
    }
    ss << "-";
    
    // flanking residue number for each residue (assume these will always be <= 9, so no delim)
    ss << MstUtils::toString(flankRes);
    name = ss.str();
    cout << "Fragment name set to:  " << name << endl;
}

vector<Structure*> seedTERM::extendMatch(seedType seed_type, const Structure* match_structure, vector<int> match_idx, vector<vector<int>> seed_segments, vector<string> seed_sec_struct, int match_number, mstreal RMSD, int match_protein_ID, bool same_res) {
    vector<Structure*> new_structures;
    
    // for each seed segment, make a new extended fragment that includes the anchor
    for (int i = 0; i < seed_segments.size(); i++) {
        
        // initialize extendedFragment structure
        Structure* extended_fragment = new Structure();
        
        // add "anchor" residues to the structure
        bool keep_chain = true;
        addResToStructure(match_idx,keep_chain,match_structure,*extended_fragment);
        
        // add "seed" residues to the structure
        keep_chain = false;
        addResToStructure(seed_segments[i],keep_chain,match_structure,*extended_fragment);
        
        // use fragment name + extendedFragment number to construct a unique name
        num_extended_fragments++;
        num_seed_residues += seed_segments[i].size();
        string ext_frag_name = extendedFragmentName(RMSD,match_protein_ID);
//        cout << "Extended fragment with name: " << ext_frag_name << endl;
        extended_fragment->setName(ext_frag_name);
        
        // add to fragmenter seeds
        Seed s;
        s.extended_fragment = extended_fragment;
        s.match_number = match_number;
        s.rmsd = RMSD;
        s.secondary_structure_classification = seed_sec_struct[i];
        s.sequence_match = same_res;
        s.num_fragment_matches = numMatches;
        
        // reassign resname to secondary structure classification
        vector<Residue*> seedChainRes = extended_fragment->getChainByID(parent->seed_chain_ID)->getResidues();
        for (int j = 0; j < seedChainRes.size(); j++) {
            Residue* R = seedChainRes[j];
            char resSSchar = s.secondary_structure_classification[j];
            string resSS = string(1,resSSchar);
            R->setName(resSS);
        }
        
        seeds.push_back(s);
        new_structures.push_back(extended_fragment);
    }
    return new_structures;
}

void seedTERM::writeExtendedFragmentstoPDB(string outDir, fstream& info_out, fstream& secstruct_out, interfaceCoverage* IC) {
    for (int i = 0; i != seeds.size(); i++) {
        Structure* ext_frag = seeds[i].extended_fragment;
        if (IC != nullptr) {
            // get seed chain atoms/residues
            Chain* seed_chain = ext_frag->getChainByID(parent->seed_chain_ID);
            
            bool overlap_peptide = IC->mapSeedToChainSubsegments(seed_chain->getAtoms(), seed_chain->getResidues(), 0, -1.0, 0);
            if (!overlap_peptide) continue;
        }
        string seed_name = outDir+ext_frag->getName()+".pdb";
        ext_frag->writePDB(seed_name);
        info_out << seed_name << endl;
    }
}

void seedTERM::writeExtendedFragmentstoBIN(fstream& info_out, fstream& secstruct_out, StructuresBinaryFile* bin, interfaceCoverage* IC) {
    for (int i = 0; i != seeds.size(); i++) {
        Structure* ext_frag = seeds[i].extended_fragment;
        
        if (IC != nullptr) {
            // get seed chain atoms/residues
            Chain* seed_chain = ext_frag->getChainByID(parent->seed_chain_ID);
            
            bool overlap_peptide = IC->mapSeedToChainSubsegments(seed_chain->getAtoms(), seed_chain->getResidues(), 0, -1.0, 0);
            if (!overlap_peptide) continue;
        }
        
        bin->appendStructure(ext_frag);
        bin->appendStructurePropertyInt("match_number",seeds[i].match_number);
        bin->appendStructurePropertyReal("match_rmsd",seeds[i].rmsd);
        bin->appendStructurePropertyReal("rmsd_adj",rmsd_adjust_factor);
        bin->appendStructurePropertyInt("seq", seeds[i].sequence_match);
        bin->appendStructurePropertyInt("num_fragment_matches", seeds[i].num_fragment_matches);
        
        //write out other information to text files
        string seed_name = ext_frag->getName();
        info_out << seed_name << endl;
        secstruct_out << seed_name << "\t" << seeds[i].secondary_structure_classification << endl;
    }
}


void seedTERM::deleteExtendedFragments() {
    for (int i = 0; i < seeds.size(); i++) {
        delete seeds[i].extended_fragment;
    }
    seeds.clear();
}

void seedTERM::searchForMatches(bool discard_match_data) {
    if (searched) return;
    else search = true;
    MstTimer timer;
    timer.start();
    
    //search for matches
    parent->F.setRMSDCutoff(RMSD_cutoff);
    parent->F.setQuery(fragmentStructure);
    parent->F.options().unsetSequenceConstraints();
    //if seq_const is true, apply
    if (parent->params.seq_const != TEParams::seqConstType::NONE) {
        if (SeqTools::isUnknown(cenRes->getName())) cout << "Warning: could not apply sequence constraint when constructing seedTERM centered around " << cenRes->getChainID() << cenRes->getNum() << " residue name: " << cenRes->getName() << " is unknown" << endl;
        Structure splitQuery = parent->F.getQuery();
        fasstSeqConstSimple seqConst(splitQuery.chainSize());
        if (parent->params.seq_const == TEParams::seqConstType::CEN_RES_ONLY) {
            const Residue& res = splitQuery.getResidue(cenResIdx);
            seqConst.addConstraint(res.getChain()->getIndex(), res.getResidueIndexInChain(), {res.getName()});
        } else if (parent->params.seq_const == TEParams::seqConstType::ALL_RES) {
            for (Residue* R : splitQuery.getResidues()) seqConst.addConstraint(R->getChain()->getIndex(), R->getResidueIndexInChain(), {R->getName()});
        } else if (parent->params.seq_const == TEParams::seqConstType::BLOSUM62) {
            if (parent->params.acceptable_aa_substitutions_path == "") MstUtils::error("seqConstType 'BLOSUM62' selected, but no path to acceptable substitions file was provided");
            for (Residue* R : splitQuery.getResidues()) {
                string singleAA = SeqTools::toSingle(R->getName());
                vector<string>& acceptable_aa_subs = parent->params.acceptable_aa_substitutions[singleAA];
                seqConst.addConstraint(R->getChain()->getIndex(), R->getResidueIndexInChain(), acceptable_aa_subs);
            }
        } else if (parent->params.seq_const == TEParams::seqConstType::HYBRID) {
            if (parent->params.acceptable_aa_substitutions_path == "") MstUtils::error("seqConstType 'HYBRID' selected, but no path to acceptable substitions file was provided");
            Residue* cenRes = &splitQuery.getResidue(cenResIdx);
            for (Residue* R : splitQuery.getResidues()) {
                if (R == cenRes) {
                    seqConst.addConstraint(cenRes->getChain()->getIndex(), cenRes->getResidueIndexInChain(), {cenRes->getName()});
                } else {
                    string singleAA = SeqTools::toSingle(R->getName());
                    vector<string>& acceptable_aa_subs = parent->params.acceptable_aa_substitutions[singleAA];
                    seqConst.addConstraint(R->getChain()->getIndex(), R->getResidueIndexInChain(), acceptable_aa_subs);
                }
            }
        }
        parent->F.options().setSequenceConstraints(seqConst);
    }
    
    // Checks for homology between the query and the match before adding to set
    tD.define(interactingRes,flankRes);
    fasstSolutionSet sols = parent->F.search();
    numMatchesPreHomologyFiltering = sols.size();
    if (parent->params.verbose) cout << "When searching found " << sols.size() << " matches. Now filtering by homology cutoff: " << parent->params.homology_cutoff << endl;
    tD.setMatches(sols,parent->params.homology_cutoff,&parent->F);
    numMatches = tD.numMatches();
    if (parent->params.verbose) cout << "After filtering have " << tD.numMatches() << " matches" << endl;
    
    timer.stop();
    construction_time = timer.getDuration();
    if (discard_match_data) tD = termData();
    if (parent->params.verbose) cout << "It took " << construction_time << " s to construct " << name << " with " << numMatches << " matches with RMSD cutoff: " << RMSD_cutoff << endl;
    searched = true;
}

int seedTERM::getSegmentNum(mstreal maxPeptideBondLen) {
    vector<int> segment_lengths = getSegmentLens(maxPeptideBondLen);
    return segment_lengths.size();
}

vector<int> seedTERM::getSegmentLens(mstreal maxPeptideBondLen) {
    vector<int> segment_lengths;
    
    vector<Residue*> residues = fragmentStructure.getResidues();
    int chain_len = 0;
    for (int i = 0; i < residues.size() - 1; i++) {
        // get last BB atom of first res and first BB atom of second
        Atom* atomC = residues[i]->findAtom("C", true);
        Atom* atomN = residues[i + 1]->findAtom("N", true);
        mstreal dist = atomC->distance(atomN);
        
        // if first residue or within the required distance add to existing segment
        if ((i == 0) || (dist < maxPeptideBondLen)) {
            chain_len += 1;
        } else {
            // otherwise, add previous segment and start new one
            segment_lengths.push_back(chain_len);
            chain_len = 1;
        }
    }
    // add final residue
    chain_len += 1;
    segment_lengths.push_back(chain_len);
    
    return segment_lengths;
}

void seedTERM::addResToStructure(vector<int> res_idx, bool keep_chain, const Structure* source_structure, Structure& recipient_structure) {
    
    for (int i = 0; i < res_idx.size(); i++) {
        // get pointer to the res in the match
        Residue* R_match = &(source_structure->getResidue(res_idx[i]));
        
        string R_frag_chain_id;
        Residue* R_new;
        if (keep_chain) {
            /* renumbers the residues to match the residue indices of the original structure
             similarly, uses the same chain name and assignments for the residues */
            
            // get pointer to the residue and the chain name in the original structure
            Residue* R_target = &(parent->getTarget()->getResidue(allResIdx[i])); //
            R_frag_chain_id = R_target->getChainID();
            
            // construct/allocate memory for new residue, based on match res
            R_new = new Residue(*R_match);
            
            // reset residue number to match the corresponding res in the original structure
            R_new->setNum(allResIdx[i]);
            
        } else {
            // as far as I know, no structures in the PDB use numbers as chain names, so I will use 0
            // set the chain ID to 0, so that residue can be added to new chain
            R_frag_chain_id = parent->seed_chain_ID;
            
            // construct/allocate memory for new residue, based on match res
            R_new = new Residue(*R_match);
            
            // reset residue number
            R_new->setNum(i);
            
        }
        
        // check chain of the residue in original protein.
        // if exists in structure, add to chain
        Chain* C = recipient_structure.getChainByID(R_frag_chain_id);
        if (C != NULL) {
            C->appendResidue(R_new);
        } else {
            // if not, initialize chain, add to structure, and add to chain
            Chain* C_new = recipient_structure.appendChain(R_frag_chain_id,false);
            C_new->appendResidue(R_new);
        }
    }
}

string seedTERM::extendedFragmentName(mstreal RMSD, int match_protein_ID) {
    stringstream ss;
    ss << name << "_";
    ss << match_protein_ID << "_";
    ss << parent->params.seed_flanking_res << "-";
    ss << RMSD << "-";
    ss << num_extended_fragments; //a unique ID
    return ss.str();
}

void TEParams::setParamsFromFile(string params_file_path) {
    // import file and split lines
    vector<string> all_lines = MstUtils::fileToArray(params_file_path);
    
    // get values from file
    for (string line : all_lines) {
        vector<string> line_split = MstUtils::split(line," ");
        if (line_split.size() != 2) MstUtils::error("Wrong number of fields ("+MstUtils::toString(line_split.size())+") in this line: "+line,"TEParams::TEParams");
        if (line_split[0] == "fragment_type") {
            string option = MstUtils::toString(line_split[1]);
            if (option == "CEN_RES") fragment_type = CEN_RES;
            else if (option == "ALL_COMBINATIONS") fragment_type = ALL_COMBINATIONS;
            else if (option == "ADAPTIVE_SIZE") fragment_type = ADAPTIVE_SIZE;
            else if (option == "ADAPTIVE_LENGTH") fragment_type = ADAPTIVE_LENGTH;
            else if (option == "ADAPTIVE_LENGTH_FIXED_RMSD") fragment_type = ADAPTIVE_LENGTH_FIXED_RMSD;
            else if (option == "COMPLEXITY_SCAN") fragment_type = COMPLEXITY_SCAN;
            else MstUtils::error("Provided value: "+option+" for option "+line_split[0]+" not recognized");
        } else if (line_split[0] == "cd_threshold") {
            cd_threshold = MstUtils::toReal(line_split[1]);
        } else if (line_split[0] == "cdSeqConst_threshold") {
            cdSeqConst_threshold = MstUtils::toReal(line_split[1]);
        } else if (line_split[0] == "int_threshold") {
            int_threshold = MstUtils::toReal(line_split[1]);
        } else if (line_split[0] == "bbInteraction_cutoff") {
            bbInteraction_cutoff = MstUtils::toReal(line_split[1]);
        } else if (line_split[0] == "max_rmsd") {
            max_rmsd = MstUtils::toReal(line_split[1]);
        } else if (line_split[0] == "flanking_res") {
            flanking_res = MstUtils::toInt(line_split[1]);
        } else if (line_split[0] == "match_req") {
            match_req = MstUtils::toInt(line_split[1]);
        } else if (line_split[0] == "adaptive_rmsd") {
            adaptive_rmsd = MstUtils::toInt(line_split[1]);
        } else if (line_split[0] == "seq_const") {
            string option = MstUtils::toString(line_split[1]);
            if (option == "NONE") seq_const = NONE;
            else if (option == "CEN_RES_ONLY") seq_const = CEN_RES_ONLY;
            else if (option == "ALL_RES") seq_const = ALL_RES;
            else if (option == "BLOSUM62") seq_const = BLOSUM62;
            else if (option == "HYBRID") seq_const = HYBRID;
            else MstUtils::error("Provided value: "+option+" for option "+line_split[0]+" not recognized");
        } else if (line_split[0] == "config_file") {
            config_file = line_split[1];
        } else if (line_split[0] == "seed_flanking_res") {
            seed_flanking_res = MstUtils::toInt(line_split[1]);
        } else if (line_split[0] == "homology_cutoff") {
            homology_cutoff = MstUtils::toReal(line_split[1]);
        } else if (line_split[0] == "allow_sidechain_clash") {
            allow_sidechain_clash = MstUtils::toInt(line_split[1]);
        } else if (line_split[0] == "freedom_cutoff") {
            freedom_cutoff = MstUtils::toReal(line_split[1]);
        } else if (line_split[0] == "relSASA_cutoff") {
            relSASA_cutoff = MstUtils::toReal(line_split[1]);
        } else if (line_split[0] == "verbose") {
            verbose = MstUtils::toInt(line_split[1]);
        } else if (line_split[0] == "acceptable_aa_substitutions_path") {
            acceptable_aa_substitutions_path = MstUtils::toString(line_split[1]);
        } else if (line_split[0] == "ignoreMissingSTRIDELabel") {
            ignoreMissingSTRIDELabel = MstUtils::toInt(line_split[1]);
        } else MstUtils::error("Option: '"+line_split[0]+"' not recognized","TEParams::TEParams");
    }
    if (acceptable_aa_substitutions_path != "") loadAASubMap();
}

void TEParams::printValues() {
    cout << "TERM Extension Parameters" << endl;
    cout << "cd_threshold " << cd_threshold << endl;
    cout << "cdSeqConst_threshold " << cdSeqConst_threshold << endl;
    cout << "int_threshold " << int_threshold << endl;
    cout << "bbInteraction_cutoff " << bbInteraction_cutoff << endl;
    cout << "max_rmsd " << max_rmsd << endl;
    cout << "flanking_res " << flanking_res << endl;
    cout << "match_req " << match_req << endl;
    cout << "adaptive_rmsd " << adaptive_rmsd << endl;
    cout << "seq_const " << seq_const << endl;
    cout << "config_file " << config_file << endl;
    cout << "seed_flanking_res " << seed_flanking_res << endl;
    cout << "homology_cutoff " << homology_cutoff << endl;
    cout << "verbose " << verbose << endl;
    cout << "ignoreMissingSTRIDELabel " << ignoreMissingSTRIDELabel << endl;
    
    if (acceptable_aa_substitutions_path != "") {
        cout << "amino acid substitution map: " << endl;
        for (auto it : acceptable_aa_substitutions) {
            cout << it.first << ": [";
            for (string aa : acceptable_aa_substitutions[it.first]) {
                cout << " " << aa;
            }
            cout << "]"<< endl;
        }
    }
}

void TEParams::loadAASubMap() {
    vector<string> lines = MstUtils::fileToArray(acceptable_aa_substitutions_path);
    for (string line : lines) {
        vector<string> split_1 = MstUtils::split(line);
        if (split_1.size() != 2) MstUtils::error("Malformed line in aaSubstition file, there should only be a single whitespace","TEParams::loadAASubMap");
        string key = split_1[0];
        vector<string> values = MstUtils::split(split_1[1],",");
        acceptable_aa_substitutions[key] = values;
    }
}


/* --------- Fragmenter --------- */
TermExtension::TermExtension(string fasstDBPath, string rotLib, vector<Residue*> selection, TEParams& _params) : sec_structure(), params(_params) {
    set_params(fasstDBPath,rotLib);
    
    cout << selection.size() << " residues selected. Setting target structure..." << endl;
    if (selection.empty()) MstUtils::error("No selected residues...","TermExtension::TermExtension");
    target_structure = selection[0]->getStructure();
    structure_name = MstSys::splitPath(target_structure->getName(),1);
    
    vCalc.setStructure(target_structure);
    
    cout << "Extracting target backbone to construct a promixity search object..." << endl;
    if (!RotamerLibrary::hasFullBackbone(*target_structure)) cout << "warning: target structure is missing backbone atoms!" << endl;
    if (params.allow_sidechain_clash) {
        target_BB_atoms = RotamerLibrary::getBackbone(*target_structure);
    } else {
        target_BB_atoms = AtomPointerVector(target_structure->getAtoms());
    }
    target_BB_structure = Structure(target_BB_atoms);
    target_PS = new ProximitySearch(target_BB_atoms, vdwRadii::maxSumRadii()/2);
    
    cout << "Building sasaCalculator..." << endl;
    sasaCalc.setStructure(target_structure);
    sasaCalc.computeSASA();
    bool relative = true;
    res2relSASA = sasaCalc.getResidueSASA(relative);
    
    cout << "Setting the binding site residues..." << endl;
    for (Residue* R : selection) {
        if (res2relSASA.at(R) > params.relSASA_cutoff) {
            bindingSiteRes.push_back(R);
        } else cout << "Residue " << *R << " with relative SASA: " << res2relSASA.at(R) << " eliminated from binding site residue set" << endl;
    }
    
    cout << "Copying atoms from structures in DB" << endl;
    //get atoms from structures in DB to minimize copying later on
    target_structures_atoms.resize(F.numTargets());
    for (int i = 0; i < F.numTargets(); i++) {
        vector<Atom*> atoms = F.getTarget(i)->getAtoms();
        vector<Atom*> atoms_copy;
        for (Atom* a: atoms) {
            Atom* a_copy = new Atom(a);
            a_copy->stripInfo();
            atoms_copy.push_back(a_copy);
        }
        target_structures_atoms[i] = atoms_copy;
    }
    
    cout << "Fragmenter construction complete." << endl;
}

void TermExtension::set_params(string fasstDBPath, string rotLib) {
    // Initialize vector to store fragments
    vector<seedTERM> all_fragments;
    
    // Import Rotamer Library
    RL_path = rotLib;
    RL.readRotamerLibrary(RL_path);
    
    // Set up FASST object for searching fragments
    foptsBase = F.options();
    F.setOptions(foptsBase);
    F.options().setRedundancyProperty("sim");
    F.options().setMinNumMatches(0); //Always allow 0 matches, min_match_number is used for the match_requirement_algorithm
    cout << "Reading FASST database... " << endl;
    MstTimer timer; timer.start();
    fasstdbPath = fasstDBPath;
    F.readDatabase(fasstdbPath,1); //only read backbone///, destroy original structure and reload as needed
    timer.stop();
    cout << "Reading the database took " << timer.getDuration() << " s... " << endl;
    
    // Check which properties it has
    cout << "Checking which properties are populated..." << endl;
    bool contProp = F.isResiduePairPropertyPopulated("cont");
    bool interferingProp = F.isResiduePairPropertyPopulated("interfering");
    bool interferedProp = F.isResiduePairPropertyPopulated("interfered");
    bool bbProp = F.isResiduePairPropertyPopulated("bb");
    
    cout << "Residue pair property: cont = " << contProp << endl;
    cout << "Residue pair property: interfering = " << interferingProp << endl;
    cout << "Residue pair property: interfered = " << interferedProp << endl;
    cout << "Residue pair property: bb = " << bbProp << endl;
    
//    if (!contProp && !interferingProp && !interferedProp && !bbProp) MstUtils::error("At least one of the following properties must be included in the DB, or else no seeds will be generated","TermExtension::set_params");
    
    // set seed/clash parameters
    minimum_seed_length = (params.seed_flanking_res*2)+1;
    extendedFragmentNumber = 0;
    
    setAAToSeqContProp();
}

TermExtension::~TermExtension() {
    delete target_PS;
    for (int i = 0; i < target_structures_atoms.size(); i++) {
        for (int j = 0; j < target_structures_atoms[i].size(); j++) delete target_structures_atoms[i][j];
    }
    for (seedTERM* f : all_fragments) delete f;
    all_fragments.clear();
}

void TermExtension::setTargetResidues(vector<Residue*> Selection) {
    bindingSiteRes = Selection;
    target_structure = bindingSiteRes[0]->getStructure();
}

void TermExtension::storeParameters(string Dir) {
    fstream params_out;
    MstUtils::openFile(params_out, Dir, fstream::out);
    params_out << "Structure name" << "\t" << target_structure->getName() << endl;
    params_out << "\tChain size" << "\t" << target_structure->chainSize() << endl;
    params_out << "\tResidue Size" << "\t" << target_structure->residueSize() << endl;
    params_out << "Databases" << endl;
    params_out << "Rotamer library" << "\t" << RL_path << endl;
    params_out << "FASST database" << "\t" << fasstdbPath << endl;
    params_out << "Fragmentation parameters" << endl;
    params_out << "\tFlanking residues" << "\t" << params.flanking_res << endl;
    params_out << "\tCD threshold" << "\t" << params.cd_threshold << endl;
    params_out << "\tINT threshold" << "\t" << params.int_threshold << endl;
    params_out << "\tBackbone interaction cutoff (Angstroms)" << "\t" << params.bbInteraction_cutoff << endl;
    params_out << "\tAdaptive RMSD" << "\t" << params.adaptive_rmsd << endl;
    params_out << "\tMatch number requirement" << "\t" << params.match_req << endl;
    params_out << "\tSequence constraint" << "\t" << params.seq_const << endl;
    params_out << "Fragment extension parameters" << endl;
    params_out << "\tSeed flanking residues" << "\t" << params.seed_flanking_res << endl;
    params_out << "\tMinimum seed length" << "\t" << minimum_seed_length << endl;
    params_out.close();
}

void TermExtension::generateFragments(bool search) {
    if (((params.fragment_type == TEParams::ADAPTIVE_SIZE)) && (!search)) {
        MstUtils::error("The algorithm used to create the selected fragment type requires that search is enabled");
    }
    
    ConFind C(&RL, *target_structure, true); // Initialize ConFind Object
    cout << MstUtils::toString(bindingSiteRes.size()) << " protein residues selected to be expanded into fragments" << endl;
    
    for (int cenResID = 0; cenResID != bindingSiteRes.size(); cenResID++) {
        Residue* cenRes = bindingSiteRes[cenResID];
        cout << "Generate fragment (" << cenResID+1 << "/" << bindingSiteRes.size() << ") centered on residue " << *(cenRes) << endl;
        
            /* CEN_RES */
        if (params.fragment_type == TEParams::CEN_RES) {
//            // Need to turn adaptive RMSD off for initial search
//            bool originalAdaptiveRMSD = params.adaptive_rmsd; // store to be reset later
//            params.adaptive_rmsd = false;
              
            // Also set the maximum number of matches to keep the memory usage down
            F.setMaxNumMatches(params.match_req);
              
            seedTERM* frag = new seedTERM(this, {cenRes}, params.flanking_res, search, params.max_rmsd);
            all_fragments.push_back(frag);
            
//            params.adaptive_rmsd = originalAdaptiveRMSD;
            F.setMaxNumMatches(-1);
            
            /* ALL_COMBINATIONS */
        } else if (params.fragment_type == TEParams::ALL_COMBINATIONS) {
            // Need to turn adaptive RMSD off for initial search
            bool originalAdaptiveRMSD = params.adaptive_rmsd; // store to be reset later
            params.adaptive_rmsd = false;
            F.setMaxNumMatches(-1);
            
            // Find the residues contacting the central residue
            vector<pair<Residue*,Residue*>> conts = generalUtilities::getContactsWith({cenRes}, C, 0, params.cd_threshold, params.int_threshold, params.bbInteraction_cutoff, params.verbose);
            vector<Residue*> contResidues = generalUtilities::getContactingResidues(conts);
            
            // Generate all possible combinations
            vector<vector<Residue*>> residueCombinations = generalUtilities::generateAllCombinationsRes(contResidues);
            
            for (vector<Residue*> resSet : residueCombinations) {
                resSet.push_back(cenRes);
                reverse(resSet.begin(),resSet.end()); //cenRes should be at front
                bool discard_data = true; //otherwise will run out of memory
                seedTERM* frag = new seedTERM(this, resSet, params.flanking_res, search, params.max_rmsd, discard_data);
                all_fragments.push_back(frag);
            }
            
            params.adaptive_rmsd = originalAdaptiveRMSD;
            
            /* ADAPTIVE_LENGTH_FIXED_RMSD */
        } else if (params.fragment_type == TEParams::ADAPTIVE_LENGTH_FIXED_RMSD) {
            if (!search) MstUtils::error("Cannot use ADAPTIVE_LENGTH_FIXED_RMSD mode to generate fragments when search is set to false","TermExtension::generateFragments");
            
            seedTERM* current_f;
            int current_flank_res = 1;
            
            // Also set the maximum number of matches to keep the memory usage down
            F.setMaxNumMatches(params.match_req);
            
            current_f = new seedTERM(this, {cenRes}, current_flank_res, search, params.max_rmsd);
            
            // If there were enough matches, try to grow by adding flanking residues
            if (params.verbose) cout << "Will now try to grow fragment by extending central segment" << endl;
            current_flank_res++;
            seedTERM* new_f;
            while (current_flank_res <= params.flanking_res) {
                // add flanking residues until params.flanking_res is reached
                new_f = new seedTERM(this, {cenRes}, current_flank_res, search, params.max_rmsd);
                
                if (new_f->getNumMatchesPreFiltering() < params.match_req) {
                    // couldn't find enough matches, store the last fragment that had sufficient matches
                    if (params.verbose) cout << "Couldn't find enough matches, store the last fragment that had sufficient matches: " << current_f->getName() << endl;
                    // delete the new fragment that didn't have enough matches
                    delete new_f;
                    
                    break;
                }
                
                if (params.verbose) cout << "Found sufficient matches to segment with " << current_flank_res << " flanking residues" << endl;
                
                // replace current fragment
                delete current_f;
                current_f = new_f;
                
                current_flank_res++;
            }
            // reset adaptive rmsd/maxNumMatches to previous state
            F.setMaxNumMatches(-1);
            
            all_fragments.push_back(current_f);
            
            /* ADAPTIVE_SIZE */
        } else if ((params.fragment_type == TEParams::ADAPTIVE_SIZE) || (params.fragment_type == TEParams::ADAPTIVE_LENGTH)) {
            if (!search) MstUtils::error("Cannot use ADAPTIVE_SIZE/ADAPTIVE_LENGTH/COMPLEXITY_SCAN mode to generate fragments when search is set to false","TermExtension::generateFragments");
            
            seedTERM* current_f;
            int current_flank_res = 1;
            
            // Need to turn adaptive RMSD off for initial search
            bool originalAdaptiveRMSD = params.adaptive_rmsd; // store to be reset later
            params.adaptive_rmsd = false;
            
            // Also set the maximum number of matches to keep the memory usage down
            F.setMaxNumMatches(params.match_req);
            
            // Search smallest reasonable fragment (three residue segment)
            if (params.verbose) cout << "Create and search smallest reasonable fragment (three residue segment)" << endl;
            current_f = new seedTERM(this, {cenRes}, current_flank_res, search, params.max_rmsd);
            
            
            
            // Check if fragment has N matches with RMSD below the adaptive cutoff
            if (params.verbose) cout << "Checking if fragment has " << params.match_req << " matches with less than " << current_f->getAdaptiveRMSD() << " adaptive RMSD" << endl;
            bool enough_matches = ((current_f->getNumMatchesPreFiltering() >= params.match_req) && (current_f->getMatch(min(params.match_req,current_f->getNumMatches())-1).getRMSD() < current_f->getAdaptiveRMSD()));
            
            if (!enough_matches) {
                // could not find enough matches to expand
                if (params.verbose) cout << "Could not find enough matches to expand, add three-residue fragment and terminate: " << current_f->getName() << endl;
                all_fragments.push_back(current_f);
                // reset adaptive rmsd to previous state
                params.adaptive_rmsd = originalAdaptiveRMSD;
                continue;
            }
            
            // Will now apply adaptive RMSD/no max num matches when searching
            params.adaptive_rmsd = true;
            
            // If there were enough matches, try to grow by adding flanking residues
            if (params.verbose) cout << "Will now try to grow fragment by extending central segment" << endl;
            current_flank_res++;
            seedTERM* new_f;
            bool should_continue = false;
            while (current_flank_res <= params.flanking_res) {
                // add flanking residues until params.flanking_res is reached
                new_f = new seedTERM(this, {cenRes}, current_flank_res, search, params.max_rmsd);
                
                if (new_f->getNumMatchesPreFiltering() < params.match_req) {
                    // couldn't find enough matches, store the last fragment that had sufficient matches
                    if (params.verbose) cout << "Couldn't find enough matches, store the last fragment that had sufficient matches: " << current_f->getName() << endl;
                    all_fragments.push_back(current_f);
                    // reset adaptive rmsd/maxNumMatches to previous state
                    params.adaptive_rmsd = originalAdaptiveRMSD;
                    F.setMaxNumMatches(-1);
                    // delete the new fragment that didn't have enough matches
                    delete new_f;
                    
                    should_continue = true;
                    break;
                }
                
                if (params.verbose) cout << "Found sufficient matches to segment with " << current_flank_res << " flanking residues" << endl;
                
                // replace current fragment
                delete current_f;
                current_f = new_f;
                
                current_flank_res++;
            }
            if (should_continue) continue;
            
            // if ADAPTIVE_LENGTH, stop before new segments are added
            if (params.fragment_type == TEParams::ADAPTIVE_LENGTH) {
                all_fragments.push_back(current_f);
                continue;
            }
            
            // The segment defined around the central segment is now the maximum possible length,
            // expand by adding new segments (residues contacting the central segment)
            
            // Turn of maxNumMatches, since we need to know how many matches to select the best fragment
            F.setMaxNumMatches(-1);
            
            // get the residues that contact the central residue
            vector<pair<Residue*,Residue*>> conts = generalUtilities::getContactsWith({cenRes}, C, 0, params.cd_threshold, params.int_threshold, params.bbInteraction_cutoff, params.verbose);
            vector<Residue*> contResidues = generalUtilities::getContactingResidues(conts);
            
            // Expansion loop: continues until there are no expansions that have N matches
            if (params.verbose) cout << "Will now try to grow fragment by adding new segments by adding contacts" << endl;
            vector<seedTERM*> expansions = generateAllExpansions(current_f,contResidues,search);
            int count = 1;
            while (!expansions.empty()) {
                if (params.verbose) cout << "Round " << count << " with " << expansions.size() << " potential expansions" << endl;
                
                // Get the best fragment
                new_f = bestFragment(expansions);
                
                // Delete all expansions that were not selected
                for (seedTERM* t : expansions) if (t != new_f) delete t;
                
                if (new_f == nullptr) {
                    // there were no fragments that satisfied the criteria, terminate
                    if (params.verbose) cout << "No viable expansions, terminate." << endl;
                    break;
                }
                
                if (params.verbose) cout << "Selecting expansion: " << new_f->getName() << endl;
                
                // Replace current fragment with new fragment
                delete current_f;
                current_f = new_f;
                
                // If a contact was added, remove it from the contResidues vector
                vector<Residue*> currentInteractingRes = current_f->getInteractingRes();
                vector<Residue*> contResiduesCopy = contResidues;
                contResidues.clear();
                for (Residue* R : contResiduesCopy) {
                    if (find(currentInteractingRes.begin(),currentInteractingRes.end(),R) == currentInteractingRes.end()) {
                        contResidues.push_back(R);
                    }
                }
                
                // Try expanding again
                expansions = generateAllExpansions(current_f,contResidues,search);
                if (params.verbose & expansions.empty()) cout << "No viable expansions, terminate." << endl;
                
                count++;
            }
            cout << "Adding: " << current_f->getName() << endl;
            all_fragments.push_back(current_f);
            
            // Reset the adaptive rmsd
            params.adaptive_rmsd = originalAdaptiveRMSD;
        }
    }
    cout << "In the end, generated the following fragments: ";
    for (seedTERM* t : all_fragments) cout << t->getName() << endl;
}


vector<Structure> TermExtension::getFragmentStructures() {
    vector<Structure> all_structures;
    for (int i = 0; i != all_fragments.size(); i++){
        all_structures.push_back(all_fragments[i]->getStructure());
    }
    return all_structures;
}

void TermExtension::writeFragmentPDBs(string outDir, string binnedDataPath) {
    oneDimBinnedData binnedData;
    if (binnedDataPath != "") {
        binnedData = oneDimBinnedData(binnedDataPath);
        binnedData.exceptOutOfRangeQueries();
    }
    
    string infoFile = outDir + "fragments.tsv";
    fstream info_out;
    MstUtils::openFile(info_out, infoFile, fstream::out);
    
    //header line
    info_out << "fragment_name\tpath\tcentral_res\tchain_id\tres_num\tnum_residue\tnum_segment\tseg_lengths\tcontacting_residues\tall_residues\tmatch_num_req\tmatch_num\tseed_num\tseed_res_num\teffective_dof\tadjusted_rmsd\trelSASA\tconstruction_time";
    if (binnedDataPath != "") info_out << "\tseed_res_per_match\tfrac_unoccupied_vol\tpredicted_seed_res_per_match\tgenerating_score";
    info_out << endl;
    
    //write line for each fragment
    for (int i = 0; i != all_fragments.size(); i++){
        seedTERM* f = all_fragments[i];
        string fragName = f->getName();
        string pdbFile = outDir + fragName + ".pdb";
        
        // write data to the file
        info_out << fragName << "\t" << pdbFile << "\t";
        info_out << f->getCenResName() << "\t" << f->getCenRes()->getChainID() << "\t";
        info_out <<  f->getCenRes()->getNum() << "\t" << f->getNumRes() << "\t" << f->getSegmentNum() << "\t";
        
        // write the segment length, for each segment
        vector<int> segment_lengths = f->getSegmentLens();
        for (int j = 0; j < segment_lengths.size(); j++) {
            info_out << segment_lengths[j];
            if (j + 1 != segment_lengths.size()) info_out << ",";
        }
        info_out << "\t";
        
        // write the contacting residues
        vector<Residue*> interacting_res = f->getInteractingRes();
        for (int j = 0; j < interacting_res.size(); j++) {
            info_out << interacting_res[j]->getChainID() << interacting_res[j]->getNum();
            if (j + 1 != interacting_res.size()) info_out << ",";
        }
        info_out << "\t";
        
        // write out all residues
        vector<int> all_res_idx = f->getResIdx();
        for (int j = 0; j < all_res_idx.size(); j++) {
            Residue* R = &target_structure->getResidue(all_res_idx[j]);
            info_out << R->getChainID() << R->getNum();
            if (j + 1 != all_res_idx.size()) info_out << ",";
        }
        info_out << "\t";
        
        // write remaining data to the file
        info_out << params.match_req << "\t" << f->getNumMatches() << "\t";
        info_out << f->getExtendedFragmentNum() << "\t" << f->getSeedResNum() << "\t";
        info_out << f->getComplexity() << "\t" << f->getSearchRMSD() << "\t";
        info_out << res2relSASA.at(f->getCenRes()) << "\t" << f->getConstructionTime();
        
        if (binnedDataPath != "") {
            info_out << "\t";
            mstreal seedResPerMatch = mstreal(f->getSeedResNum())/mstreal(f->getNumMatches()+1);
            mstreal fracUnoccupiedVol = vCalc.fracUnoccupiedVolume(f->getCenRes());
            mstreal expectedCount = binnedData.getVal(fracUnoccupiedVol);
            mstreal generating_score = -log2((seedResPerMatch+0.1)/(expectedCount+0.1));
            info_out << seedResPerMatch << "\t" << fracUnoccupiedVol << "\t";
            info_out << expectedCount << "\t" << generating_score;
            for (Atom* A : f->getCenRes()->getAtoms()) A->setB(generating_score);
        }
        
        info_out << endl;
        
        // Now store the structure itself
        f->getStructure().writePDB(pdbFile);
        
    }
    info_out.close();
    
    //write the target structure, with the B-factor set to the seed residue generating score, if it was calculated
    target_structure->writePDB(outDir+structure_name+"seedResidueGeneratingScore.pdb");
}

void TermExtension::writeFragmentClassification(string outDir) {
    string infoFile = outDir + "fragment_res_classification.tsv";
    fstream info_out;
    MstUtils::openFile(info_out, infoFile, fstream::out);
    
    info_out << "fragment_name\tchain_id\tres_num\tclassification" << endl;
        
    for (seedTERM* f : all_fragments) {
        vector<int> all_res_idx = f->getResIdx();
        for (int res_idx : all_res_idx){
            Residue* R = &target_structure->getResidue(res_idx);
            string classification = sec_structure.classifyResidue(R);
            info_out << f->getName() << "\t" << R->getChainID() << "\t" << R->getNum() << "\t" << classification << endl;
        }
    }
}


void TermExtension::extendFragmentsandWriteStructures(seedTERM::seedType option, string outDir) {
    //  if ((file_type != "pdb") && (file_type != "bin")) MstUtils::error("File type not recognized: " +file_type);
    
    MstTimer timer;
    timer.start();
    
    string matchContactsFile = outDir + "matchContacts.tsv";
    string binFile = outDir + "extendedfragments.bin";
    string infoFile = outDir + "extendedfragments.info";
    string secStrucFile = outDir + "seed_secondary_structure.info";
    fstream match_out, info_out, secstruct_out;
    
    //open files to write
    if (params.verbose) MstUtils::openFile(match_out, matchContactsFile, fstream::out);
    MstUtils::openFile(info_out, infoFile, fstream::out);
    MstUtils::openFile(secstruct_out, secStrucFile, fstream::out);
    StructuresBinaryFile* bin = new StructuresBinaryFile(binFile,false); //open in write mode, start new file
    
    //write header
    if (params.verbose) match_out << "targetID\tcenResID\tmatchRes\tcontactRes\tseqConstContactRes\tinterferingRes\tinterferedRes\tbbInteractionRes\tselectedRes\tclashFilteredRes\tseedNames" << endl;
    
    cout << "Extending " << all_fragments.size() << " fragments total, each with no more than " << params.match_req << " matches" << endl;
    
    // Now, iterate over each target that has at least one match
    for (int frag_id = 0; frag_id < all_fragments.size(); frag_id++) {
        seedTERM* frag = all_fragments[frag_id];
        int num_new_structures = extendFragment(frag, option, bin, match_out,  info_out, secstruct_out);
        extendedFragmentNumber += num_new_structures;
    }
    if (params.verbose) match_out.close();
    info_out.close();
    secstruct_out.close();
    
    timer.stop();
    cout << "Seed creation took " << timer.getDuration() << " s" << endl;
    delete bin;
}

int TermExtension::extendFragment(seedTERM* frag, seedTERM::seedType option, StructuresBinaryFile* bin, fstream& match, fstream& info, fstream& sec_struct) {
    int num_extendedfragments = 0;
    fasstSolutionSet& matches = frag->getMatches();
    vector<Structure*> fragment_extensions;
    int N_matches;
    if (params.match_req < 0) N_matches = matches.size();
    else N_matches = min(matches.size(),params.match_req);
    cout << "Extending fragment: " << frag->getName() << endl;
    cout << "matches size: " << matches.size() << " match req: " << params.match_req << endl;
    //matches are sorted by RMSD, so take the top N_matches
    for (int sol_id = 0; sol_id < N_matches; sol_id++) {
        //    if (params.verbose) cout << "Fragment: " << frag->getName() << " Solution ID: " << sol_id << endl;
        fasstSolution& sol = matches[sol_id];
        
        //transform the match for this specific solution
        int target_idx = sol.getTargetIndex();
        if (params.verbose) match << target_idx << "\t";
        Structure* transformed_match_structure = F.getTarget(target_idx);
        sol.getTransform().apply(transformed_match_structure);
        
        //// Look for contacts to the central residue in the match protein and create seeds by the
        //// specified method
        
        vector<int> match_idx = F.getMatchResidueIndices(sol,FASST::matchType::REGION);
        
        // get the idx of the residues that 1) contact the central residue of the match and 2) could
        // potentially be used to construct a seed (flanking residues don't overlap). Return these +
        // the union of their flanking residues to be used to construct a seed.
        
        string R_query_aaName = frag->getCenRes()->getName();

        
        vector<int> seed_res_idx = identifySeedResidueIdx(frag, transformed_match_structure, match_idx, sol.getTargetIndex(), R_query_aaName, match);
        
        
        // for each residue, check for potential clashes between its backbone atoms and the target
        // proteins backbone atoms
        vector<int> filtered_seed_res_idx = getNonClashingResidueIdx(seed_res_idx,transformed_match_structure,match);
        
        // separate all contiguous sequences of residues (where there is no gap) that are greater than
        // the minimum_seed_length into their own vector
        vector<vector<int>> seed_segments = getSeedSegments(filtered_seed_res_idx, transformed_match_structure);
        
        vector<string> seed_sec_structure = classifySegmentsUsingSTRIDELabels(target_idx,seed_segments);
        
        //check if central residue of the match has same sequence as query
        string R_match_aaName = transformed_match_structure->getResidue(match_idx[frag->getCenResIdx()]).getName();
        bool same_res = (R_match_aaName == R_query_aaName);
        
        //use the match_idx/extension_idx to make a fragment extension
        vector<Structure*> new_structures = frag->extendMatch(option,transformed_match_structure,match_idx,seed_segments,seed_sec_structure,sol_id,sol.getRMSD(),target_idx,same_res);
        
        for (Structure* extfrag : new_structures) {
            match << extfrag->getName() << "," << endl;
        }
        
        if (storeAll) {
            transformed_match_structure->writePDB(wholeMatchProteinDir+frag->getName()+"_"+MstUtils::toString(target_idx)+"_wholeMatch.pdb");
            for (Structure* S : new_structures) {
                S->writePDB(extFragDir+S->getName()+"_extendedFragment.pdb");
                Structure seed = Structure(*S->getChainByID(seed_chain_ID));
                seed.writePDB(seedDir+S->getName()+"_seed.pdb");
            }
        }
        
        frag->writeExtendedFragmentstoBIN(info,sec_struct,bin,IC);
        
        num_extendedfragments += new_structures.size();
        
        frag->deleteExtendedFragments();
        
        //reset the coordinates of the target structure
        generalUtilities::copyAtomCoordinates(transformed_match_structure, target_structures_atoms[target_idx]);
        
        if (params.verbose) match << endl;
    }
    cout << "Fragment " << frag->getName() << " generated " << frag->getExtendedFragmentNum() << " extended fragments!" << endl;
    return num_extendedfragments;
}

vector<seedTERM*> TermExtension::generateAllExpansions(seedTERM* fragToExpand, vector<Residue*> contactingRes, bool search) {
    vector<seedTERM*> fragExpansions;
    
    // Find new contacts to the central residue and add to fragment as new segment
    vector<Residue*> interactingRes = fragToExpand->getInteractingRes();
    vector<Residue*> tryAdding =  MstUtils::setdiff(contactingRes,interactingRes);
    
    for (Residue* newRes : tryAdding) {
        vector<Residue*> newInteractingRes = interactingRes;
        newInteractingRes.insert(newInteractingRes.end(),newRes);
        
         // Construct new fragment and add if not already equivalent to a member
        seedTERM* new_f = new seedTERM(this,newInteractingRes,params.flanking_res,false,params.max_rmsd);
        fragExpansions.push_back(new_f);;
    }
    return fragExpansions;
}

seedTERM* TermExtension::bestFragment(vector<seedTERM*> fragExpansions) {
    // The goal is to eliminate uneccessary searches. The fragments are sorted by size
    
    if (fragExpansions.empty()) return nullptr;
    sort(fragExpansions.begin(),fragExpansions.end());
    int currentSize = -1;
    vector<seedTERM*> currentGroup;
    for (auto it = fragExpansions.rbegin(); it != fragExpansions.rend(); it++) {
        // check if this fragment is smaller than the last one and we have some viable fragments in
        // current group
        if (((*it)->getNumRes() <= currentSize) && (!currentGroup.empty())) {
            // select the fragment with the most matches and terminate
            sort(currentGroup.begin(),currentGroup.end());
            return currentGroup.back();
        }
        // search for matches, and if has enough, add to currentGroup
        (*it)->searchForMatches();
        if ((*it)->getNumMatchesPreFiltering() >= params.match_req) currentGroup.push_back(*it);
        currentSize = (*it)->getNumRes();
    }
    if (!currentGroup.empty()) {
        // select the fragment with the most matches and terminate
        sort(currentGroup.begin(),currentGroup.end());
        return currentGroup.back();
    }
    return nullptr;
}


vector<int> TermExtension::identifySeedResidueIdx(const seedTERM* frag, const Structure* match_structure, vector<int> match_idx, int fasst_target_index, string& cenResAA, fstream& match) {
    
    // Identify the id of the central residue of the fragment, in the match
    int cenResIdx = frag->getCenResIdx();
    int cenResMatchIdx = match_idx[cenResIdx];
    
    if (params.verbose) match << cenResMatchIdx << "\t";
    if (params.verbose) match << joinString(match_idx, ",") << "\t";

    
    // Identify residues interacting with central res in match protein using EACH contact type
    map<int, mstreal> all_contacts = F.getResiduePairProperties(fasst_target_index,"cont",cenResMatchIdx);
    
    // Sequence constrained contacts are special, and each is stored as a unique property
    map<int, mstreal> all_contacts_seqconst;
    if (aaToSeqContProp.count(cenResAA) == 1) {
        all_contacts_seqconst = F.getResiduePairProperties(fasst_target_index,aaToSeqContProp[cenResAA],cenResMatchIdx);
    }
    
    map<int, mstreal> all_interfered = F.getResiduePairProperties(fasst_target_index,"interfered",cenResMatchIdx);
    map<int, mstreal> all_interfering = F.getResiduePairProperties(fasst_target_index,"interfering",cenResMatchIdx);
    map<int, mstreal> all_bbInteraction =  F.getResiduePairProperties(fasst_target_index,"bb",cenResMatchIdx);
    
    // Filter each contact type by its respective value, add to filtered contacts if new
    //CD
    vector<int> filtered_conts;
    vector<int> contacts_for_type;
    for (auto it = all_contacts.begin(); it != all_contacts.end(); it++) {
        if (it->second > params.cd_threshold) {
            filtered_conts.push_back(it->first);
            if (params.verbose) contacts_for_type.push_back(it->first);
        }
    }
    if (params.verbose) {
        match << joinString(contacts_for_type,",") << "\t";
        contacts_for_type.clear();
    }
    
    //CD seq const
    for (auto it = all_contacts_seqconst.begin(); it != all_contacts_seqconst.end(); it++) {
        if ((it->second > params.cdSeqConst_threshold) && (find(filtered_conts.begin(),filtered_conts.end(),it->first) == filtered_conts.end())) {
            filtered_conts.push_back(it->first);
            if (params.verbose) contacts_for_type.push_back(it->first);
        }
    }
    if (params.verbose) {
        match << joinString(contacts_for_type,",") << "\t";
        contacts_for_type.clear();
    }
    //Int (residues that are interfered with by central res)
    for (auto it = all_interfered.begin(); it != all_interfered.end(); it++) {
        if ((it->second > params.int_threshold) && (find(filtered_conts.begin(),filtered_conts.end(),it->first) == filtered_conts.end())) {
            filtered_conts.push_back(it->first);
            if (params.verbose) contacts_for_type.push_back(it->first);
        }
    }
    if (params.verbose) {
        match << joinString(contacts_for_type,",") << "\t";
        contacts_for_type.clear();
    }
    //Int (residues that interfere with the central res)
    for (auto it = all_interfering.begin(); it != all_interfering.end(); it++) {
        if ((it->second > params.int_threshold) && (find(filtered_conts.begin(),filtered_conts.end(),it->first) == filtered_conts.end())) {
            filtered_conts.push_back(it->first);
            if (params.verbose) contacts_for_type.push_back(it->first);
        }
    }
    if (params.verbose) {
        match << joinString(contacts_for_type,",") << "\t";
        contacts_for_type.clear();
    }
    //bb
    for (auto it = all_bbInteraction.begin(); it != all_bbInteraction.end(); it++) {
        if ((it->second < params.bbInteraction_cutoff) && (find(filtered_conts.begin(),filtered_conts.end(),it->first) == filtered_conts.end())) {
            filtered_conts.push_back(it->first);
            if (params.verbose) contacts_for_type.push_back(it->first);
        }
    }
    if (params.verbose) {
        match << joinString(contacts_for_type,",") << "\t";
        contacts_for_type.clear();
    }
    
    // Search for contacts that do not already appear in the residue indices of the match (or their flank)
    // if slow, replace this
    vector<int> seed_conts;
    for (int i = 0; i < filtered_conts.size(); i++) {
        int cont_res = filtered_conts[i];
        // note the extra +/- 1. We want to ensure that these residues can act as "central residues"
        // for the case where neither the target/peptide is being extended. So we need to ensure that
        // it's potential flanking residues do not overlap with the match residues AND that they
        // are not covalently bonded to match residues.
        if ((find(match_idx.begin(), match_idx.end(), cont_res+params.flanking_res+1) == match_idx.end()) && (find(match_idx.begin(), match_idx.end(), cont_res-params.flanking_res-1) == match_idx.end())) {
            seed_conts.push_back(cont_res);
        }
    }
    
    //  if (params.verbose) cout << "identifySeedResidueIdx() seed_conts residue number " << seed_conts.size() << endl;
    
    // Construct a vector that contains these contacts + the union of their flanking residues
    // sort the vector
    sort(seed_conts.begin(),seed_conts.end());
    
    // if slow, replace this
    vector<int> seed_res;
    for (int i = 0; i < seed_conts.size(); i++) {
        // In the case that a contact is near the boundary - avoid creating a residue index for a
        // residue that does not exist
        
        // Note: this can technically give flanking residues that are on another chain, because this
        // index is residue level. This is accounted for after clash filtering, when the actual seed
        // is constructed
        int min_idx = max(seed_conts[i]-params.seed_flanking_res, 0);
        int max_idx = min(seed_conts[i]+params.seed_flanking_res, match_structure->residueSize() - 1);
        for (int j = min_idx; j <= max_idx; j++) {
            // only add if not already added
            if (find(seed_res.begin(),seed_res.end(),j) == seed_res.end()) seed_res.push_back(j);
        }
    }
    if (params.verbose) match << joinString(seed_res,",") << "\t";
    return seed_res;
}

vector<int> TermExtension::getNonClashingResidueIdx(vector<int> seed_res_idx, const Structure* match_structure, fstream& match) {
    vector<int> filtered_res_idx;
    for (int i = 0; i < seed_res_idx.size(); i++) {
        Residue* R = &match_structure->getResidue(seed_res_idx[i]);
        AtomPointerVector R_bb_atoms = AtomPointerVector(RotamerLibrary::getBackbone(R));
        bool clash = isClash(*target_PS, target_BB_atoms, R_bb_atoms);
        if (!clash) filtered_res_idx.push_back(seed_res_idx[i]);
    }
    if (params.verbose) {
        for (int filt_res : filtered_res_idx) match << filt_res << (filt_res != filtered_res_idx.back() ? "," : "");
        match << "\t";
            };
    return filtered_res_idx;
}

vector<vector<int>> TermExtension::getSeedSegments(vector<int> seed_res_idx, const Structure* match_structure) {
    vector<vector<int>> all_segments;
    /* As the residues are sorted by index, all that needs to be checked is 1) whether
     the previous index is within the flanking region of this residue. If false, a new segment
     is required. First residue is always added to the segment.
     
     If segment is completed, but falls below the minimum length requirement, it will not be added*/
    int previous_index;
    Chain* previous_chain;
    vector<int> segment;
    for (int i = 0; i < seed_res_idx.size(); i++) {
        int index = seed_res_idx[i];
        Chain* chain = match_structure->getResidue(index).getChain();
        if ((i == 0) || ((index-1 == previous_index) && (chain == previous_chain))) {
            segment.push_back(index);
        } else {
            // the previous segment is now completed, add if longer than minimum seed length
            if (segment.size() >= minimum_seed_length) all_segments.push_back(segment);
            // remove the residues in the segment
            segment.clear();
            // add the current residue to the fresh segment
            segment.push_back(index);
        }
        previous_index = index;
        previous_chain = chain;
    }
    // there will always be a remaining segment to add
    if (segment.size() >= minimum_seed_length) all_segments.push_back(segment);
    return all_segments;
}

vector<string> TermExtension::classifySegmentsInMatchProtein(Structure* S, vector<vector<int>> seed_segments) {
    vector<string> classifications;
    for (vector<int> segment : seed_segments) classifications.push_back(sec_structure.classifyResInStruct(S, segment));
    return classifications;
}

vector<string> TermExtension::classifySegmentsUsingSTRIDELabels(int target_idx, vector<vector<int>> seed_segments) {
    vector<string> classifications;
    for (vector<int> segment : seed_segments) {
        string segment_ss = "";
        string res_stride;
        for (int res_idx : segment) {
            res_stride = F.getResidueStringProperty(target_idx, "stride", res_idx);
            if ((res_stride.empty())||(res_stride == " ")) {
                if (params.ignoreMissingSTRIDELabel) res_stride = "C";
                else MstUtils::error("Residue STRIDE label missing for residue in DB","TermExtension::classifySegmentsUsingSTRIDELabels");
            }
            segment_ss += res_stride;
        }
        classifications.push_back(segment_ss);
    }
    return classifications;
}

void TermExtension::setAAToSeqContProp() {
    cout << "Set the AAToSeqContProp map" << endl;
    vector<string> aaNames = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
        "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL", "ALA"}; // all but GLY and PRO
    for (string aa : aaNames) {
        aaToSeqContProp[aa] = "cont"+aa;
        cout << "key: " << aa << ", value: " << "cont"+aa << endl;
    }
}

//vector<Structure*> TermExtension::getExtendedFragments() {
//  vector<Structure*> all_seeds;
//  for (int i = 0; i < all_fragments.size(); i ++) {
//    vector<Structure*> seeds = all_fragments[i]->getExtendedFragments();
//    all_seeds.insert(all_seeds.end(), seeds.begin(), seeds.end());
//  }
//  return all_seeds;
//}

void nearbySeedScore::scoreResidues() {
    for (Residue* R : selectedRes) selRes2NearbySeedResCount[R] = 0.0;

    /**
     To save on memory, load seeds in batches and for each seed residue, count nearby protein residues. Since fragments may have different
     numbers of matches, adjust by a proportionality factor so that the counts are comparable.
     */
    while (seeds.hasNext()) {
        vector<Structure*> seed_batch = seeds.next();
        for (Structure* seed : seed_batch) {
            mstreal numMatches = seeds.getStructurePropertyInt("num_fragment_matches", seed->getName());
            for (Residue* R : seed->getResidues()) {
                CartesianPoint CA_coor = R->findAtom("CA")->getCoor();
                vector<int> sel_CA_nearby = psCA.getPointsWithin(CA_coor, 0, cR);
                for (int point : sel_CA_nearby) {
                    Atom* A = selCAatoms[point];
                    // adjust depending on the number of matches
                    mstreal proportion = maxNumMatches/numMatches;
                    selRes2NearbySeedResCount[A->getResidue()] += 1.0*proportion;
                }
            }
        }
    }

    // Divide by the total number of matches and find score
    for (auto const& it : selRes2NearbySeedResCount) {
        Residue* R = it.first;
        mstreal val = it.second;
        selRes2NearbySeedResCount[R] = val/maxNumMatches;
        selRes2FracUnoccupiedVol[R] = vCalc.fracUnoccupiedVolume(R);
        mstreal expectedCount = 0.0;
        if (binnedDataPath != "") expectedCount = binnedData.getVal(selRes2FracUnoccupiedVol[R]);
        selRes2Score[R] = -log2((selRes2NearbySeedResCount[R]+0.1)/(expectedCount+0.1));
    }
}

void nearbySeedScore::writeResidueInfo(string prefix, set<Residue*> bindingSiteRes) {
    fstream info_out;
    string path = prefix + "nearbySeedResidueScore.csv";
    MstUtils::openFile(info_out,path,fstream::out);
    
    // set all B-factors to unreasonably high number
    for (Residue* R : target->getResidues()) for (Atom* A : R->getAtoms()) A->setB(100.0);
    
    info_out << "chain_id,res_num,aa,frac_unoccupied_volume,nearby_seed_res,predicted_nearby_seed_res,hotspot_score,binding_site_res" << endl;
    
    for (Residue* R : selectedRes) {
        mstreal expectedCount = 0.0;
        if (binnedDataPath != "") expectedCount = binnedData.getVal(selRes2FracUnoccupiedVol[R]);
        
        info_out << R->getChainID() << "," << R->getNum() << ",";
        info_out << R->getName() << ",";
        info_out << selRes2FracUnoccupiedVol[R] << ",";
        info_out << selRes2NearbySeedResCount[R] << ",";
        info_out << expectedCount << ",";
        info_out << selRes2Score[R] << ",";
        info_out << ((bindingSiteRes.find(R) != bindingSiteRes.end()) ? 1 : 0) << endl;
        
        // Set residue atoms B-factor to score
        for (Atom* A : R->getAtoms()) A->setB(selRes2Score[R]);
    }
    info_out.close();
    
    string pdbName = prefix + "nearbySeedResidueScore.pdb";
    target->writePDB(pdbName);
}

