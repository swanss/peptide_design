//
//  generateFragments.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 2/1/19.
//

#include "termextension.h"

using namespace MST;

/* --------- Fragmenter --------- */
TermExtension::TermExtension(string fasstDBPath, string rotLib, vector<Residue*> selection, bool _verbose) : sec_structure() {
    verbose = _verbose;
    set_params(fasstDBPath,rotLib);
    
    if (verbose) cout << selection.size() << " residues selected. Setting target structure..." << endl;
    if (selection.empty()) MstUtils::error("No selected residues...");
    target_structure = selection[0]->getStructure();
    structure_name = MstSys::splitPath(target_structure->getName(),1);
    
    if (verbose) cout << "Extracting target backbone to construct a promixity search object..." << endl;
    if (!RotamerLibrary::hasFullBackbone(*target_structure)) cout << "warning: target structure is missing backbone atoms!" << endl;
    target_BB_atoms = RotamerLibrary::getBackbone(*target_structure);
    target_BB_structure = Structure(target_BB_atoms);
    target_PS = new ProximitySearch(target_BB_atoms, vdwRadii::maxSumRadii()/2);
    
    if (verbose) cout << "Setting the binding site residues..." << endl;
    bindingSiteRes = selection;
    
    if (verbose) cout << "Copying structures from DB" << endl;
    //get structures from DB to minimize copying
    for (int i = 0; i < F.numTargets(); i++) {
        target_structures.push_back(new Structure(*F.getTarget(i)));
    }
    
    cout << "Fragmenter construction complete." << endl;
}


void TermExtension::set_params(string fasstDBPath, string rotLib) {
    // Store the structure/parameters
    cd_threshold = .01;
    int_threshold = .01;
    bbInteraction_cutoff = 3.5;
    max_rmsd = 1.2;
    flanking_res = 2;
    match_req = -1; //only considered if value is > 0
    adaptive_rmsd = true;
    seq_const = false;
    
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
    if (verbose) cout << "Reading FASST database... ";
    MstTimer timer; timer.start();
    fasstdbPath = fasstDBPath;
    F.readDatabase(fasstdbPath,1); //only read backbone
    timer.stop();
    if (verbose) cout << "Reading the database took " << timer.getDuration() << " s... " << endl;
    
    // set seed/clash parameters
    seed_flanking_res = 2;
    minimum_seed_length = (seed_flanking_res*2)+1;
    extendedFragmentNumber = 0;
    
}

TermExtension::~TermExtension() {
    delete target_PS;
    for (int i = 0; i < target_structures.size(); i++) delete target_structures[i];
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
    params_out << "\tFlanking residues" << "\t" << flanking_res << endl;
    params_out << "\tCD threshold" << "\t" << cd_threshold << endl;
    params_out << "\tINT threshold" << "\t" << int_threshold << endl;
    params_out << "\tBackbone interaction cutoff (Angstroms)" << "\t" << bbInteraction_cutoff << endl;
    params_out << "\tAdaptive RMSD" << "\t" << adaptive_rmsd << endl;
    if (match_req != 0) params_out << "\tMatch number requirement" << "\t" << match_req << endl;
    params_out << "Fragment extension parameters" << endl;
    params_out << "\tSeed flanking residues" << "\t" << seed_flanking_res << endl;
    params_out << "\tMinimum seed length" << "\t" << minimum_seed_length << endl;
    params_out.close();
}

void TermExtension::generateFragments(fragType option, bool search) {
    if (((option == MATCH_NUM_REQ_SIZE) || (option == MATCH_NUM_REQ_CUTOFF) || (option == COMPLEXITY_SCAN)) && (!search)) {
        MstUtils::error("The algorithm used to create the selected fragment type requires that search is enabled");
    }
    
    ConFind C(&RL, *target_structure, true); // Initialize ConFind Object
    cout << MstUtils::toString(bindingSiteRes.size()) << " protein residues selected to be expanded into fragments" << endl;
    
    for (int cenResID = 0; cenResID != bindingSiteRes.size(); cenResID++) {
        Residue* cenRes = bindingSiteRes[cenResID];
        if (verbose) cout << "Generate fragment (" << cenResID+1 << "/" << bindingSiteRes.size() << ") centered on residue " << *(cenRes) << endl;
        
            /* CEN_RES */
        if (option == CEN_RES) {
            seedTERM* frag = new seedTERM(this, {cenRes}, search, max_rmsd, flanking_res);
            all_fragments.push_back(frag);
            
            /* CONTACT */
        } else if (option == CONTACT) {
            vector<pair<Residue*,Residue*>> conts = generalUtilities::getContactsWith({cenRes}, C, 0, cd_threshold, int_threshold, bbInteraction_cutoff, verbose);
            vector<Residue*> allRes = generalUtilities::getContactingResidues(conts);
            // place center residue at the front of the vector of all residues
            allRes.insert(allRes.begin(), cenRes);
            seedTERM* frag = new seedTERM(this, allRes, search, max_rmsd, flanking_res);
            all_fragments.push_back(frag);
            
            /* ALL_COMBINATIONS */
        } else if (option == ALL_COMBINATIONS) {
            // get the residues that contact the central residue
            vector<pair<Residue*,Residue*>> conts = generalUtilities::getContactsWith({cenRes}, C, 0, cd_threshold, int_threshold, bbInteraction_cutoff, verbose);
            vector<Residue*> contResidues = generalUtilities::getContactingResidues(conts);
            
            //generate all possible combinations of this set of residues
            vector<vector<Residue*>> contacting_res_combinations = generalUtilities::generateAllCombinationsRes(contResidues);
            
            if (verbose) cout << "There are " << contacting_res_combinations.size() << " possible combinations of residues interacting with the center residue" << endl;
            
            for (int combination = 0; combination < contacting_res_combinations.size(); combination++) {
                vector<Residue*>& contacting_res = contacting_res_combinations[combination];
                vector<Residue*> all_res;
                all_res.reserve(contacting_res.size() + 1);
                all_res.insert(all_res.end(), bindingSiteRes[cenResID]);
                all_res.insert(all_res.end(), contacting_res.begin(), contacting_res.end());
                seedTERM* frag = new seedTERM(this, all_res, search, max_rmsd, flanking_res);
                all_fragments.push_back(frag);
            }
            /* MATCH_GUIDED (find proper RMSD cutoff) */
        } else if (option == MATCH_NUM_REQ_CUTOFF) {
            // Begin by making a fragment with the binding site residue alone
            seedTERM* self_f = new seedTERM(this, {cenRes}, search, max_rmsd, flanking_res);
            if (verbose) cout << "Constructed fragment with default cutoff " << max_rmsd << endl;
            
            //if this has the required number of matches, stop
            if (self_f->getNumMatches() > match_req) {
                if (verbose) cout << "Fragment has sufficient matches, add and continue" << endl;
                all_fragments.push_back(self_f);
            } else {
                //if there are insufficient matches, raise the cutoff
                if (verbose) cout << "Fragment has insufficient matches, try raising cutoff" << endl;
                mstreal new_rmsd_cutoff = max_rmsd;
                int new_flanking_res = flanking_res;
                seedTERM* self_f = nullptr;
                while (new_flanking_res > 0) {
                    while (new_rmsd_cutoff <= 1.0) {
                        cout << "Try raising the cutoff to " << new_rmsd_cutoff << endl;
                        delete self_f;
                        self_f = new seedTERM(this, {cenRes}, search, new_rmsd_cutoff, new_flanking_res);
                        if (self_f->getNumMatches() > match_req) {
                            if (verbose) cout << "Fragment has sufficient matches, add and continue" << endl;
                            all_fragments.push_back(self_f);
                            goto end;
                        }
                        new_rmsd_cutoff += 0.1;
                    }
                    //decrement the flanking residue number, making a shorter segment
                    cout << "Try decreasing number of flanking residues to " << new_flanking_res-1 << " (and reset cutoff)" << endl;
                    new_flanking_res--;
                    new_rmsd_cutoff = max_rmsd;
                }
                if (verbose) cout << "Was not able to find sufficient matches to fragment within allowed RMSD cutoff range, adding anyways" << endl;
                all_fragments.push_back(self_f);
                end:;
            }
            /* MATCH_GUIDED (find proper fragment size) */
        } else if (option == MATCH_NUM_REQ_SIZE || option == COMPLEXITY_SCAN) {
            /*
             MATCH_NUM_REQ_SIZE: Try adding each residue from the contacting residues to the central one
             and select the largest fragment (highest number of residues) that 1) has at least the required
             number of matches and 2) more total residues than the previous fragment. If there is a tie,
             take the fragment with more matches. If nothing can be added, store the current fragment.
             
             COMPLEXITY_SCAN: Similar to above, but after each cycle of selecting a residue to add,
             store the fragment, also, continue until there are no interacting residues remaining
             Try combining each contacting residue with the central one to generate a fragment and see if there are
             enough matches. If so, create a fragment to be grown later.
             */
            
            // Begin by making a fragment with the binding site residue alone
            seedTERM* self_f = new seedTERM(this, {cenRes}, search, max_rmsd, flanking_res);
            
            // get the residues that contact the central residue
            vector<pair<Residue*,Residue*>> conts = generalUtilities::getContactsWith({cenRes}, C, 0, cd_threshold, int_threshold, bbInteraction_cutoff, verbose);
            vector<Residue*> contResidues = generalUtilities::getContactingResidues(conts);
            
            // if no contacts - then just add the fragment and continue
            if (contResidues.empty() || option == COMPLEXITY_SCAN) {
                all_fragments.push_back(self_f);
                if (contResidues.empty()) continue;
            }
            
            //Begin loop (I realize this is a nightmare to read through)
            vector<Residue*> remConts = contResidues;
            seedTERM* current_fragment = self_f;
            while (!remConts.empty()) {
                map<Residue*,seedTERM*> expansions;
                vector<Residue*> expansion_res;
                
                // Try adding each of the remaining residues to the current residue
                for (int i = 0; i < remConts.size(); i++) {
                    Residue* R = remConts[i];
                    if (verbose) cout << "Try adding contact " << *R << "..." << endl;
                    vector<Residue*> new_res = current_fragment->getInteractingRes();
                    new_res.push_back(R);
                    seedTERM* new_f = new seedTERM(this, new_res, search, max_rmsd, flanking_res);
                    
                    // If fragment has at least the minimum number of matches AND more res than current, keep
                    if (new_f->getNumMatches() >= match_req && new_f->getNumRes() > current_fragment->getNumRes()) {
                        if (verbose) cout << new_f->getName() << " added to the potential expansions" << endl;
                        expansions[R] = new_f;
                        expansion_res.push_back(R);
                    } else {
                        delete new_f;
                    }
                }
                
                // Now decide whether to stop OR select an expansion to continue expanding in the next iteration
                if (expansions.empty()) {
                    // No contacts can be added, break loop
                    if (verbose) cout << "No potential expansions... process complete" << endl;
                    // Add current fragment in match_num_req mode, as this is the final fragment
                    if (option == MATCH_NUM_REQ_SIZE) all_fragments.push_back(current_fragment);
                    if ((option == MATCH_NUM_REQ_SIZE) && (verbose)) cout << "In the end, " << current_fragment->getName() << " will be added (MATCH_NUM_REQ_SIZE)" << endl;
                    // In the COMPLEXITY_SCAN mode, this should already have been added last cycle
                    break;
                } else if (expansion_res.size() == 1) {
                    // Only one contact passing the criteria, add
                    if (option == COMPLEXITY_SCAN) all_fragments.push_back(expansions[expansion_res.front()]);
                    if (option == COMPLEXITY_SCAN && verbose) cout << expansions[expansion_res.front()]->getName() << " added to fragments (COMPLEXITY_SCAN)" << endl;
                    // Remove the added contact
                    remConts = MstUtils::setdiff(remConts,{expansion_res.front()});
                    // If this is was the last of the remaining contacts, loop will end after this iteration, so this is the final fragment
                    if (remConts.empty() && (option == MATCH_NUM_REQ_SIZE)) all_fragments.push_back(expansions[expansion_res.front()]);
                    if (remConts.empty() && (option == MATCH_NUM_REQ_SIZE) && verbose) cout << "In the end, " << expansions[expansion_res.front()]->getName() << " will be added (MATCH_NUM_REQ_SIZE)" << endl;
                    // If not final iteration, set current fragment ()
                    if (option == MATCH_NUM_REQ_SIZE) delete current_fragment; //only the final fragment should be kept
                    if (!remConts.empty()) current_fragment = expansions[expansion_res.front()];
                    if (!remConts.empty() && verbose) cout << expansions[expansion_res.front()]->getName() << " is now the current fragment" << endl;
                    /* Don't break because--while unlikely--if there are remaining contacts, it's possible
                     that adding another contact in the context of the existing structure could somehow
                     *increase* number of matches. This is due to the influence of the adaptive RMSD cutoff,
                     which increases with the complexity of the fragment */
                } else if (expansion_res.size() > 1) {
                    // Of the available fragments, find the largest and add
                    sort(expansion_res.begin(), expansion_res.end(), [&expansions](Residue* i, Residue* j) {return expansions[i]->getNumRes() > expansions[j]->getNumRes(); });
                    if (expansions[expansion_res[0]]->getNumRes() > expansions[expansion_res[1]]->getNumRes()) { //should avoid a segfault by checking vector size earlier...
                        // The first fragment is the largest
                        if (option == COMPLEXITY_SCAN) all_fragments.push_back(expansions[expansion_res[0]]);
                        if (option == COMPLEXITY_SCAN && verbose) cout << expansions[expansion_res[0]]->getName() << " added to fragments as it is the largest (COMPLEXITY_SCAN)" << endl;
                        // Remove the added contact
                        remConts = MstUtils::setdiff(remConts,{expansion_res[0]});
                        // set current fragment
                        if (option == MATCH_NUM_REQ_SIZE) delete current_fragment; //only the final fragment should be kept
                        current_fragment = expansions[expansion_res[0]];
                        if (verbose) cout << expansions[expansion_res[0]]->getName() << " is now the current fragment" << endl;
                        //delete all other fragments that were not chosen
                        for (int j = 1; j < expansion_res.size(); j++) {
                            delete expansions[expansion_res[j]];
                        }
                    } else {
                        // If the first fragment is not the largest, find all equal size fragments and select the one with the most matches
                        if (verbose) cout << "Multiple fragments have the same residue number, breaking tie by comparing match number" << endl;
                        vector<Residue*> equal_size_expansion_res;
                        equal_size_expansion_res.push_back(expansion_res.front()); //compare the remaining to the first expansion
                        for (int j = 1; j < expansion_res.size(); j++) {
                            if (expansions[expansion_res[j]]->getNumRes() == expansions[expansion_res.front()]->getNumRes()) {
                                equal_size_expansion_res.push_back(expansion_res[j]);
                            } else {
                                delete expansions[expansion_res[j]]; //delete all other fragments that were not chosen
                            }
                        }
                        sort(equal_size_expansion_res.begin(), equal_size_expansion_res.end(), [&expansions](Residue* i, Residue* j) {return expansions[i]->getNumMatches() > expansions[j]->getNumMatches(); });
                        
                        if (option == COMPLEXITY_SCAN) all_fragments.push_back(expansions[equal_size_expansion_res[0]]);
                        if (option == COMPLEXITY_SCAN && verbose) cout << expansions[equal_size_expansion_res[0]]->getName() << " added to fragments as it has the most matches (COMPLEXITY_SCAN)" << endl;
                        
                        // Remove the added contact
                        remConts = MstUtils::setdiff(remConts,{equal_size_expansion_res[0]});
                        
                        // set current fragment
                        if (option == MATCH_NUM_REQ_SIZE) delete current_fragment; //only the final fragment should be kept
                        current_fragment = expansions[equal_size_expansion_res[0]];
                        if (verbose) cout << expansions[equal_size_expansion_res[0]]->getName() << " is now the current fragment" << endl;
                        //delete all other fragments that were not chosen
                        for (int j = 1; j < equal_size_expansion_res.size(); j++) {
                            delete expansions[equal_size_expansion_res[j]];
                        }
                    }
                }
            }
        }
    }
}


vector<Structure> TermExtension::getFragmentStructures() {
    vector<Structure> all_structures;
    for (int i = 0; i != all_fragments.size(); i++){
        all_structures.push_back(all_fragments[i]->getStructure());
    }
    return all_structures;
}

void TermExtension::writeFragmentPDBs(string outDir) {
    //begin by writing the target structure (in case it was renamed/renumbered)
    target_structure->writePDB(outDir+structure_name+".pdb");
    string infoFile = outDir + "fragments.tsv";
    fstream info_out;
    MstUtils::openFile(info_out, infoFile, fstream::out);
    
    //header line
    info_out << "fragment_name\tpath\tcentral_res\tnum_residue\tnum_segment\tseg_lengths\tcontacting_residues\tall_residues\tmatch_num_req\tmatch_num\tseed_num\teffective_dof\tadjusted_rmsd\tconstruction_time" << endl;
    
    //write line for each fragment
    for (int i = 0; i != all_fragments.size(); i++){
        seedTERM* f = all_fragments[i];
        string fragName = f->getName();
        string pdbFile = outDir + fragName + ".pdb";
        
        // write data to the file
        info_out << fragName << "\t" << pdbFile << "\t" << f->getCenResName() << "\t" << f->getNumRes() << "\t" << f->getSegmentNum() << "\t";
        
        // write the segment length, for each segment
        vector<int> segment_lengths = f->getSegmentLens();
        for (int j = 0; j < segment_lengths.size(); j++) {
            info_out << segment_lengths[j];
            if (j + 1 != segment_lengths.size()) info_out << ",";
        }
        info_out << "\t";
        
        // write the contacting residues
        vector<Residue*> interacting_res = f->getInteractingRes();
        for (int j = 1; j < interacting_res.size(); j++) {
            info_out << interacting_res[j]->getChainID() << interacting_res[j]->getNum();
            if (j + 1 != interacting_res.size()) info_out << ",";
        }
        info_out << "\t";
        
        // write out all residues
        vector<int> all_res_idx = f->getResIdx();
        for (int j = 1; j < all_res_idx.size(); j++) {
            Residue* R = &target_structure->getResidue(all_res_idx[j]);
            info_out << R->getChainID() << R->getNum();
            if (j + 1 != all_res_idx.size()) info_out << ",";
        }
        info_out << "\t";
        
        // write remaining data to the file
        info_out << match_req << "\t" << f->getNumMatches() << "\t" << f->getExtendedFragmentNum() << "\t" << f->getComplexity() << "\t" << f->getAdjRMSD() << "\t" << f->getConstructionTime() << endl;
        
        // Now store the structure itself
        f->getStructure().writePDB(pdbFile);
    }
    info_out.close();
}

void TermExtension::writeFragmentClassification(string outDir) {
    string infoFile = outDir + "fragment_res_classification.tsv";
    fstream info_out;
    MstUtils::openFile(info_out, infoFile, fstream::out);
    
    info_out << "fragment_name\tchain_id\tres_num\tclassification" << endl;
    
    secondaryStructureClassifier classifier;
    
    for (seedTERM* f : all_fragments) {
        vector<int> all_res_idx = f->getResIdx();
        for (int res_idx : all_res_idx){
            Residue* R = &target_structure->getResidue(res_idx);
            string classification = classifier.classifyResidue(R);
            info_out << f->getName() << "\t" << R->getChainID() << "\t" << R->getNum() << "\t" << classification << endl;
        }
    }
}


void TermExtension::extendFragmentsandWriteStructures(seedTERM::seedType option, string outDir) {
    //  if ((file_type != "pdb") && (file_type != "bin")) MstUtils::error("File type not recognized: " +file_type);
    
    MstTimer timer;
    timer.start();
    
    string binFile = outDir + "extendedfragments.bin";
    string infoFile = outDir + "extendedfragments.info";
    string secStrucFile = outDir + "seed_secondary_structure.info";
    fstream info_out, secstruct_out;
    
    //open files to write
    MstUtils::openFile(info_out, infoFile, fstream::out);
    MstUtils::openFile(secstruct_out, secStrucFile, fstream::out);
    StructuresBinaryFile* bin = new StructuresBinaryFile(binFile,false); //open in write mode, start new file
    
    cout << "Extending " << all_fragments.size() << " fragments total" << endl;
    
    // Now, iterate over each target that has at least one match
    for (int frag_id = 0; frag_id < all_fragments.size(); frag_id++) {
        seedTERM* frag = all_fragments[frag_id];
        int num_new_structures = extendFragment(frag, option, bin, info_out, secstruct_out);
        extendedFragmentNumber += num_new_structures;
    }
    info_out.close();
    secstruct_out.close();
    
    timer.stop();
    cout << "Seed creation took " << timer.getDuration() << " s" << endl;
    delete bin;
}

int TermExtension::extendFragment(seedTERM* frag, seedTERM::seedType option, StructuresBinaryFile* bin, fstream& info, fstream& sec_struct) {
    int num_extendedfragments = 0;
    fasstSolutionSet& matches = frag->getMatches();
    vector<Structure*> fragment_extensions;
    int N_matches;
    if (match_req <= 0) N_matches = matches.size();
    else N_matches = min(matches.size(),match_req);
    //matches are sorted by RMSD, so take the top N_matches
    for (int sol_id = 0; sol_id < N_matches; sol_id++) {
        //    if (verbose) cout << "Fragment: " << frag->getName() << " Solution ID: " << sol_id << endl;
        fasstSolution& sol = matches[sol_id];
        
        //transform the match for this specific solution
        int target_idx = sol.getTargetIndex();
        Structure* transformed_match_structure = target_structures[target_idx];
        sol.getTransform().apply(transformed_match_structure);
        
        //// Look for contacts to the central residue in the match protein and create seeds by the
        //// specified method
        
        vector<int> match_idx = F.getMatchResidueIndices(sol,FASST::matchType::REGION);
        
        // get the idx of the residues that 1) contact the central residue of the match and 2) could
        // potentially be used to construct a seed (flanking residues don't overlap). Return these +
        // the union of their flanking residues to be used to construct a seed.
        vector<int> seed_res_idx = identifySeedResidueIdx(frag, transformed_match_structure, match_idx, sol.getTargetIndex());
        
        
        // for each residue, check for potential clashes between its backbone atoms and the target
        // proteins backbone atoms
        vector<int> filtered_seed_res_idx = getNonClashingResidueIdx(seed_res_idx,transformed_match_structure);
        
        // separate all contiguous sequences of residues (where there is no gap) that are greater than
        // the minimum_seed_length into their own vector
        vector<vector<int>> seed_segments = getSeedSegments(filtered_seed_res_idx, transformed_match_structure);
        
        vector<string> seed_sec_structure = classifySegmentsInMatchProtein(transformed_match_structure,seed_segments);
        
        //check if central residue of the match has same sequence as query
        Residue* R_match = &transformed_match_structure->getResidue(match_idx[frag->getCenResIdx()]);
        Residue* R_query = frag->getCenRes();
        bool same_res = (R_match->getName() == R_query->getName());
        
        //use the match_idx/extension_idx to make a fragment extension
        vector<Structure*> new_structures = frag->extendMatch(option,transformed_match_structure,match_idx,seed_segments,seed_sec_structure,sol_id,sol.getRMSD(),same_res);
        
        frag->writeExtendedFragmentstoBIN(info,sec_struct,bin);
        
        num_extendedfragments += new_structures.size();
        
        frag->deleteExtendedFragments();
        
        //reset the coordinates of the target structure
        generalUtilities::copyAtomCoordinates(transformed_match_structure, F.getTarget(target_idx));
    }
    if (verbose) cout << "Fragment " << frag->getName() << " generated " << frag->getExtendedFragmentNum() << " extended fragments!" << endl;
    return num_extendedfragments;
}

vector<int> TermExtension::identifySeedResidueIdx(const seedTERM* frag, const Structure* match_structure, vector<int> match_idx, int fasst_target_index) {
    
    // Identify the id of the central residue of the fragment, in the match
    int cenResIdx = frag->getCenResIdx();
    int cenResMatchIdx = match_idx[cenResIdx];
    
    // Identify residues interacting with central res in match protein using EACH contact type
    map<int, mstreal> all_contacts = F.getResiduePairProperties(fasst_target_index,"cont",cenResMatchIdx);
    map<int, mstreal> all_interfered = F.getResiduePairProperties(fasst_target_index,"interfered",cenResMatchIdx);
    map<int, mstreal> all_interfering = F.getResiduePairProperties(fasst_target_index,"interfering",cenResMatchIdx);
    map<int, mstreal> all_bbInteraction =  F.getResiduePairProperties(fasst_target_index,"bb",cenResMatchIdx);
    
    // Filter each contact type by its respective value, add to filtered contacts if new
    //CD
    vector<int> filtered_conts;
    for (auto it = all_contacts.begin(); it != all_contacts.end(); it++) {
        if (it->second > cd_threshold) filtered_conts.push_back(it->first);
    }
    //Int (residues that are interfered with by central res)
    for (auto it = all_interfered.begin(); it != all_interfered.end(); it++) {
        if ((it->second > int_threshold) && (find(filtered_conts.begin(),filtered_conts.end(),it->first) == filtered_conts.end())) filtered_conts.push_back(it->first);
    }
    //Int (residues that interfere with the central res)
    for (auto it = all_interfering.begin(); it != all_interfering.end(); it++) {
        if ((it->second > int_threshold) && (find(filtered_conts.begin(),filtered_conts.end(),it->first) == filtered_conts.end())) filtered_conts.push_back(it->first);
    }
    //bb
    for (auto it = all_bbInteraction.begin(); it != all_bbInteraction.end(); it++) {
        if ((it->second < bbInteraction_cutoff) && (find(filtered_conts.begin(),filtered_conts.end(),it->first) == filtered_conts.end())) filtered_conts.push_back(it->first);
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
        if ((find(match_idx.begin(), match_idx.end(), cont_res+flanking_res+1) == match_idx.end()) && (find(match_idx.begin(), match_idx.end(), cont_res-flanking_res-1) == match_idx.end())) {
            seed_conts.push_back(cont_res);
        }
    }
    
    //  if (verbose) cout << "identifySeedResidueIdx() seed_conts residue number " << seed_conts.size() << endl;
    
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
        int min_idx = max(seed_conts[i]-seed_flanking_res, 0);
        int max_idx = min(seed_conts[i]+seed_flanking_res, match_structure->residueSize() - 1);
        for (int j = min_idx; j <= max_idx; j++) {
            // only add if not already added
            if (find(seed_res.begin(),seed_res.end(),j) == seed_res.end()) seed_res.push_back(j);
        }
    }
    return seed_res;
}

vector<int> TermExtension::getNonClashingResidueIdx(vector<int> seed_res_idx, const Structure* match_structure) {
    vector<int> filtered_res_idx;
    for (int i = 0; i < seed_res_idx.size(); i++) {
        Residue* R = &match_structure->getResidue(seed_res_idx[i]);
        AtomPointerVector R_bb_atoms = AtomPointerVector(RotamerLibrary::getBackbone(R));
        bool clash = isClash(*target_PS, target_BB_atoms, R_bb_atoms);
        if (!clash) filtered_res_idx.push_back(seed_res_idx[i]);
    }
//      if (verbose) cout << seed_res_idx.size() - filtered_res_idx.size()  << " residues removed due to clashes" << endl;
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

//vector<Structure*> TermExtension::getExtendedFragments() {
//  vector<Structure*> all_seeds;
//  for (int i = 0; i < all_fragments.size(); i ++) {
//    vector<Structure*> seeds = all_fragments[i]->getExtendedFragments();
//    all_seeds.insert(all_seeds.end(), seeds.begin(), seeds.end());
//  }
//  return all_seeds;
//}

/* --------- Fragment --------- */
seedTERM::seedTERM() {
    // compiler error thrown unless I include this
}

seedTERM::seedTERM(TermExtension* FragmenterObj, vector<Residue*> allRes, bool _search, mstreal max_rmsd, int flank) : RMSD_cutoff(max_rmsd), flank(flank) {
    MstTimer timer;
    timer.start();
    
    // set up the fragment
    parent = FragmenterObj;
    cenRes = allRes[0];
    cenResName = cenRes->getChainID() + MstUtils::toString(cenRes->getNum());
    interactingRes = allRes;
    num_extended_fragments = 0;
    
    //generate the structure
    cenResIdx = TERMUtils::selectTERM(interactingRes, fragmentStructure, flank, &(allResIdx))[0];
    
    //get the complexity
    effective_degrees_of_freedom = generalUtilities::fragDegreesOfFreedom(allResIdx, *parent->getTarget(), 15);
    rmsd_adjust_factor = 1.0/sqrt(fragmentStructure.atomSize()/effective_degrees_of_freedom);
    
    //generate a name that is unique to the fragment
    //NOTE: this name uses the residue numbers from the parent structure
    setName();
    
    //reset the residue numbers to the residue indices of the parent structure, so that they can be
    //easily added to a fusion topology.
    for (int i = 0; i < fragmentStructure.residueSize(); i++) {
        Residue* R = &fragmentStructure.getResidue(i);
        R->setNum(allResIdx[i]);
        
        num_extended_fragments = 0;
    }
    
    search = _search;
    if (search) {
        //search for matches
        FragmenterObj->F.setQuery(fragmentStructure);
        if (FragmenterObj->getAdaptiveRMSD()) RMSD_cutoff = RMSDCalculator::rmsdCutoff(allResIdx, *parent->getTarget(), RMSD_cutoff, 15); // NOTE: the residues in the same chain of the original structure are not by default assumed to be independent
        FragmenterObj->F.setRMSDCutoff(RMSD_cutoff);
        
        FragmenterObj->F.options().unsetSequenceConstraints();
        //if seq_const is true, apply
        if (parent->seq_const) {
            Structure splitQuery = FragmenterObj->F.getQuery();
            fasstSeqConstSimple seqConst(splitQuery.chainSize());
            const Residue& res = fragmentStructure.getResidue(cenResIdx);
            seqConst.addConstraint(res.getChain()->getIndex(), res.getResidueIndexInChain(), {res.getName()});
            FragmenterObj->F.options().setSequenceConstraints(seqConst);
        }
        
        matches = FragmenterObj->F.search();
    }
    
    timer.stop();
    construction_time = timer.getDuration();
    cout << "It took " << construction_time << " s to construct " << name << " with " << matches.size() << " matches with RMSD cutoff: " << RMSD_cutoff << endl;
};

seedTERM::seedTERM(const seedTERM& F) {
    parent = F.parent;
    fragmentStructure = F.fragmentStructure;
    cenRes = F.cenRes;
    cenResName = F.cenResName;
    interactingRes = F.interactingRes;
    allResIdx = F.allResIdx;
    cenResIdx = F.cenResIdx; //index of the index of the central residue
    flank = F.flank;
    matches = F.matches;
    search = F.search;
    seeds = F.seeds;
    RMSD_cutoff = F.RMSD_cutoff;
    effective_degrees_of_freedom = F.effective_degrees_of_freedom;
    rmsd_adjust_factor = F.rmsd_adjust_factor;
    name = F.name;
    construction_time = F.construction_time;
    num_extended_fragments = F.num_extended_fragments;
}

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
    
    // protein flanking residue number
    ss << MstUtils::toString(flank);
    name = ss.str();
    cout << "Fragment name set to:  " << name << endl;
}

vector<Structure*> seedTERM::extendMatch(seedType seed_type, const Structure* match_structure, vector<int> match_idx, vector<vector<int>> seed_segments, vector<string> seed_sec_struct, int match_number, mstreal RMSD, bool same_res) {
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
        string ext_frag_name = extendedFragmentName(RMSD);
//        cout << "Extended fragment with name: " << ext_frag_name << endl;
        extended_fragment->setName(ext_frag_name);
        
        // add to fragmenter seeds
        Seed s;
        s.extended_fragment = extended_fragment;
        s.match_number = match_number;
        s.rmsd = RMSD;
        s.secondary_structure_classification = seed_sec_struct[i];
        s.sequence_match = same_res;
        seeds.push_back(s);
        new_structures.push_back(extended_fragment);
    }
    return new_structures;
}

void seedTERM::writeExtendedFragmentstoPDB(string outDir, fstream& info_out, fstream& secstruct_out) {
    for (int i = 0; i != seeds.size(); i++) {
        Structure* ext_frag = seeds[i].extended_fragment;
        string seed_name = outDir+ext_frag->getName()+".pdb";
        ext_frag->writePDB(seed_name);
        info_out << seed_name << endl;
    }
}

void seedTERM::writeExtendedFragmentstoBIN(fstream& info_out, fstream& secstruct_out, StructuresBinaryFile* bin) {
    for (int i = 0; i != seeds.size(); i++) {
        Structure* ext_frag = seeds[i].extended_fragment;
        bin->appendStructure(ext_frag);
        bin->appendStructurePropertyInt("match_number",seeds[i].match_number);
        bin->appendStructurePropertyReal("match_rmsd",seeds[i].rmsd);
        bin->appendStructurePropertyReal("rmsd_adj",rmsd_adjust_factor);
        bin->appendStructurePropertyInt("seq", seeds[i].sequence_match);
        
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
            R_frag_chain_id = "0";
            
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

string seedTERM::extendedFragmentName(mstreal RMSD) {
    stringstream ss;
    ss << name << "_";
    ss << parent->seed_flanking_res << "-";
    ss << RMSD << "-";
    ss << num_extended_fragments; //a unique ID
    return ss.str();
}

