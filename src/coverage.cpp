#include "coverage.h"

using namespace MST;

/* --------- sortedBins --------- */
void sortedBins::insert(seedSubstructureInfo& entry, mstreal val) {
    if (val > max_val) return;
    if (bins.empty()) MstUtils::error("Tried to insert an entry to an improperly constructed instance of sortedBins");
    bins[val2Bin(val)].push_back(entry);
    sorted = false;
    empty = false;
}

vector<seedSubstructureInfo> sortedBins::getBinByValue(mstreal val) {
    if (val > max_val) {
        vector<seedSubstructureInfo> empty;
        return empty;
    }
    return bins[val2Bin(val)];
};

vector<seedSubstructureInfo> sortedBins::getLowestValuePopulatedBin() {
    for (auto& bin: bins) {
        if (!bin.empty()) return bin;
    }
    return bins.front();
}

seedSubstructureInfo sortedBins::getLowestRMSDSeed(set<string> seedNames) {
    for (auto& bin : bins){
        if (!bin.empty()) {
            if (seedNames.count(bin.front().structure_name) == 0) return bin.front();
            else continue;
        }
    }
    return seedSubstructureInfo();
}

vector<tuple<mstreal,mstreal,long>> sortedBins::getSeedsByBin() {
    vector<tuple<mstreal,mstreal,long>> output;
    mstreal min_boundary = 0;
    for (int i = 0; i < bins.size(); i++) {
        output.push_back(tuple<mstreal,mstreal,long> (min_boundary,min_boundary+interval,bins[i].size()));
        min_boundary += interval;
    }
    return output;
}

//vector<tuple<mstreal,mstreal,long>> sortedBins::getNumFragmentsCoveringContact(int R_prot_idx) {
//  vector<tuple<mstreal,mstreal,long>> output;
//  mstreal min_boundary = 0;
//  for (int i = 0; i < bins.size(); i++) {
//    //find all seedSubstructures that have a protein residue matching the query
//    long count = 0;
//    for (int j = 0; j < bins[i].size(); j++) {
//      if (bins[i][j].protein_res_idx.count(R_prot_idx) > 0) count++;
//    }
//    output.push_back(tuple<mstreal,mstreal,long> (min_boundary,min_boundary+interval,count));
//
//    min_boundary += interval;
//  }
//  return output;
//}

vector<seedSubstructureInfo> sortedBins::getAllSeeds() {
    vector<seedSubstructureInfo> all_seeds;
    for (vector<seedSubstructureInfo> seed_bin : bins) {
        all_seeds.insert(all_seeds.end(),seed_bin.begin(),seed_bin.end());
    }
    return all_seeds;
}

void sortedBins::sortBins() {
    for (int i = 0; i < bins.size(); i++) {
        sort(bins[i].begin(),bins[i].end());
    }
    sorted = true;
}

void sortedBins::reset() {
    for (int i = 0; i < bins.size(); i++) {
        bins[i].clear();
    }
}

void sortedBins::buildBins() {
    
    int bin_number = ceil(max_val/interval);
    
    bins.resize(bin_number);
    
}

///* --------- allChainSegments --------- */
//allChainSegments::allChainSegments(Chain* peptide, Structure* _target, int max_segment_length, mstreal _max_rmsd, string _s_cid, RotamerLibrary* RL) : target(_target), max_rmsd(_max_rmsd), s_cid(_s_cid), rotLib(RL) {
//  int peptide_length = peptide->residueSize();
//  vector<Atom*> peptide_atoms = peptide->getAtoms();
//  max_allowable_segment_length = min(peptide_length,max_segment_length);
//
//  //resize the outer vector based on the number of allowable segment lengths
//  chainSubsegments.resize(max_allowable_segment_length);
//  chainSubsegmentsBins.resize(max_allowable_segment_length);
//
//  for (int segment_length = 1; segment_length <= max_allowable_segment_length; segment_length++) {
//    int total_num_segments = peptide_length - segment_length + 1;
//
//    //resize the inner vector based on the number of allowable subsegments
//    chainSubsegments[segment_length-1].resize(total_num_segments);
//    chainSubsegmentsBins[segment_length-1].resize(total_num_segments);
//
//    for (int segment_position = 0; segment_position < total_num_segments; segment_position++) {
//      int start_position = segment_position*4;
//      int end_position = (segment_length+segment_position)*4;
//      vector<Atom*> peptide_segment(peptide_atoms.begin()+start_position,peptide_atoms.begin()+end_position);
//
//      sortedBins bin(max_rmsd);
//
//      chainSubsegments[segment_length-1][segment_position] = peptide_segment;
//      chainSubsegmentsBins[segment_length-1][segment_position] = bin;
//    }
//  }

//  RL_path = rotLibPath;
//  rotLib.readRotamerLibrary(RL_path);
//  cout << "rotamer library is loaded? " << rotLib.isLoaded() << endl;
//}



/* --------- interfaceCoverage --------- */
interfaceCoverage::interfaceCoverage(Structure& S, string p_cid, string _RL_path) : binFilePath("") {
    // intial construction
    setParams(_RL_path);
    
    //extract backbone
    if (!RotamerLibrary::hasFullBackbone(S)) {
        MstUtils::error("Structure is missing backbone atoms");
    }
  
    original_complex = &S;
    complex = Structure(RotamerLibrary::getBackbone(S));
    complex.setName(S.getName());
    
    peptide_chain = complex.getChainByID(p_cid);
    peptide_first_residue_index_in_structure = peptide_chain->getResidue(0).getResidueIndex();
    
    //define coverage elements
    cout << "Identifying structural elements..." << endl;
    peptide_residues = peptide_chain->getResidues();
    cout << "Total number of peptide residue elements: " << peptide_residues.size() << endl;
    
    all_types_contacts = defineContacts(complex, peptide_residues);
    
    contact_residues = generalUtilities::mergeContactLists({all_types_contacts["contact"],all_types_contacts["interfering"],all_types_contacts["interfered"],all_types_contacts["bbinteraction"]});
    
    cout << "Total number of contact residue elements: " << contact_residues.size() << endl;
    
    for (pair<Residue*,Residue*> cont : contact_residues) cout << cont.first->getChainID() << cont.first->getNum() << "-" << cont.second->getChainID() << cont.second->getNum() << " ";
    cout << endl;
    
    //construct target (protein chains only) and get binding site residues
    cout << "Constructing target structure and identifying binding site residues" << endl;
    prepareForTERMExtension();
    
    //define and construct all subsegments of the peptide chain
    cout << "Constructing peptide subsegments" << endl;
    int peptide_length = peptide_chain->residueSize();
    vector<Residue*> peptide_residues = peptide_chain->getResidues();
    vector<Atom*> peptide_atoms = peptide_chain->getAtoms();
    max_allowable_segment_length = min(peptide_length,max_seed_length);
    
    //resize the outer vector based on the number of allowable segment lengths
    chainSubsegments.resize(max_allowable_segment_length);
    chainSubsegmentsBins.resize(max_allowable_segment_length);
    chainResidueSubsegments.resize(max_allowable_segment_length);
    
    for (int segment_length = 1; segment_length <= max_allowable_segment_length; segment_length++) {
        int total_num_segments = peptide_length - segment_length + 1;
        
        //resize the inner vector based on the number of allowable subsegments
        chainSubsegments[segment_length-1].resize(total_num_segments);
        chainSubsegmentsBins[segment_length-1].resize(total_num_segments);
        chainResidueSubsegments[segment_length-1].resize(total_num_segments);
        
        for (int segment_position = 0; segment_position < total_num_segments; segment_position++) {
            // allocator is [start,end)
            int start_position = segment_position*4;
            int end_position = (segment_length+segment_position)*4;
            vector<Residue*> peptide_residue_segment(peptide_residues.begin()+segment_position,peptide_residues.begin()+segment_position+segment_length);
            vector<Atom*> peptide_segment(peptide_atoms.begin()+start_position,peptide_atoms.begin()+end_position);
            
            sortedBins bin(max_rmsd);
            
            chainResidueSubsegments[segment_length-1][segment_position] = peptide_residue_segment;
            chainSubsegments[segment_length-1][segment_position] = peptide_segment;
            chainSubsegmentsBins[segment_length-1][segment_position] = bin;
        }
    }
}

interfaceCoverage::~interfaceCoverage() {
    if (binFilePath != "") {
        delete seeds;
    }
    delete target;
}

void interfaceCoverage::setSeeds(string _binFilePath) {
    if (seeds != nullptr) {
        resetBins();
    }
    if (binFilePath != "") {
        delete seeds;
    }
    binFilePath = _binFilePath;
    seeds = new StructuresBinaryFile(binFilePath);
    seeds->scanFilePositions();
    seeds->reset();
}

void interfaceCoverage::setSeeds(StructuresBinaryFile* bin) {
    if (seeds != nullptr) {
        resetBins();
    }
    if (binFilePath != "") {
        delete seeds;
    }
    seeds = bin;
    seeds->scanFilePositions();
}

void interfaceCoverage::findCoveringSeeds() {
    resetBins();
    seeds->reset();
    while (seeds->hasNext()) {
        Structure* extended_fragment = seeds->next();
        
        //check if structure meets criteria
        int match_number = seeds->getStructurePropertyInt("match_number",extended_fragment->getName());
        if ((match_number > match_number_cutoff) || (match_number_cutoff == 0)) {
            delete extended_fragment;
            continue;
        }
        
        bool sequence_match = seeds->getStructurePropertyInt("seq",extended_fragment->getName());
        if (seq_const && !sequence_match) {
            //this match had a different amino acid at the central position, ignore
            delete extended_fragment;
            continue;
        }
        mstreal match_rmsd = seeds->getStructurePropertyReal("match_rmsd",extended_fragment->getName());
        if (match_rmsd_cutoff != 0.0) {
            mstreal rmsd_adjust_factor = seeds->getStructurePropertyReal("rmsd_adj",extended_fragment->getName());
            if (match_rmsd >= match_rmsd_cutoff*rmsd_adjust_factor) {
                //rmsd_adjust_factor is always in the range (0,1)
                //this match had an rmsd (after adjustment) that didn't pass the cutoff, ignore
                delete extended_fragment;
                continue;
            }
        }
        
        Chain* seed_C = extended_fragment->getChainByID(seed_chain_id);
        if (seed_C == NULL) MstUtils::error("Structures in binary file are missing seed chains: "+seed_chain_id,"interfaceCoverage::findCoveringSeeds()");
        mapSeedToChainSubsegments(getBackboneAtoms(seed_C),seed_C->getResidues(),match_number,match_rmsd,sequence_match);
        delete extended_fragment;
    }
}

bool interfaceCoverage::mapSeedToChainSubsegments(vector<Atom*> seed_atoms, vector<Residue*> seed_residues, int match_number, mstreal match_rmsd, bool sequence_match, bool only_check_if_aligned) {
    bool aligned = false;
    int seed_length = seed_residues.size();
    int max_allowable_seed_segment_length = min(seed_length,max_allowable_segment_length);
    
    //try all (allowable) segment lengths
    for (int segment_length = 1; segment_length <= max_allowable_seed_segment_length; segment_length++){
        int total_num_segments = seed_length - segment_length + 1;
        
        //slide across the seed
        for (int seed_position = 0; seed_position < total_num_segments; seed_position++) {
            int start_position = seed_position*4;
            int end_position = (seed_position+segment_length)*4;
            vector<Atom*> seed_segment(seed_atoms.begin()+start_position,seed_atoms.begin()+end_position);
            
            //slide across the peptide
            for (int peptide_position = 0; peptide_position < chainSubsegments[segment_length-1].size(); peptide_position++) {
                vector<Atom*>& peptide_segment = chainSubsegments[segment_length-1][peptide_position];
                
                mstreal rmsd = rmsd_calc.rmsd(peptide_segment,seed_segment);
                
                //check if passes the cutoff
                if (rmsd < getMaxRMSDForSegLength(max_rmsd,segment_length)) {
                    aligned = true;
                    if (!only_check_if_aligned) {
                        Atom* A = seed_segment[0];
                        string structure_name = A->getStructure()->getName();
                        string chain_ID = A->getChain()->getID();
                        
                        //compute the mean cosine angle between residue normal vectors
                        vector<Residue*> seed_residue_segment(seed_residues.begin()+seed_position,seed_residues.begin()+seed_position+segment_length);
                        mstreal alignment_cos_angle = generalUtilities::avgCosAngleBetweenSegments(seed_residue_segment,chainResidueSubsegments[segment_length-1][peptide_position]);
                        
                        seedSubstructureInfo info(structure_name,chain_ID,seed_position,segment_length,match_number,match_rmsd,sequence_match,rmsd,alignment_cos_angle);
                        
                        chainSubsegmentsBins[segment_length-1][peptide_position].insert(info, rmsd);
                    }
                }
            }
        }
    }
    return aligned;
}

void interfaceCoverage::writePeptideResidues(string outDir) {
    //residues
    cout << "residue info..." << endl;
    fstream out;
    string output_path = outDir + "residues.tsv";
    MstUtils::openFile(out, output_path, fstream::out);
    
    //write the header
    out << "peptide_position\tresidue_number\tchain\tb_factor\tnum_contacting_prot_res" << endl;
    
    vector<Residue*> peptide_residues_in_contact;
    for (pair<Residue*,Residue*> cont : contact_residues) peptide_residues_in_contact.push_back(cont.first);
    
    for (Residue* R : peptide_residues) {
        Atom* Ca = R->findAtom("CA");
        int num_contacting_prot_res = count(peptide_residues_in_contact.begin(),peptide_residues_in_contact.end(),R);
        out << R->getResidueIndexInChain() << "\t";
        out << R->getNum() << "\t";
        out << R->getChainID() << "\t";
        out << Ca->getB() << "\t";
        out << num_contacting_prot_res << endl;
    }
    out.close();
}

void interfaceCoverage::writeContacts(string outDir, set<pair<Residue*,Residue*>> provided_contact_residues, map<string,contactList> provided_all_types_contacts) {
    
    // decide whether to use provided contacts or those already defined in the object
    set<pair<Residue*,Residue*>> _contact_residues;
    map<string,contactList> _all_types_contacts;
    if (provided_contact_residues.empty()) {
        _contact_residues = contact_residues;
        _all_types_contacts = all_types_contacts;
        cout << "No contacts provided, using predefined set: " << _contact_residues.size() << endl;
    } else {
        _contact_residues = provided_contact_residues;
        if (provided_all_types_contacts.empty()) MstUtils::error("If contacting residues are passed as arguments, a map defining all types of contacts must be as well","interfaceCoverage::writeContacts");
        _all_types_contacts = provided_all_types_contacts;
         cout << "Contacts provided: " << _contact_residues.size() << endl;
    }
    
    cout << "contact info..." << endl;
    fstream out;
    string output_path = outDir + "contacts.tsv";
    MstUtils::openFile(out, output_path, fstream::out);
    
    out << "peptide_residue_number\tpeptide_chain\tprotein_residue_number\tprotein_chain\tcontact_name\tcontact\tinterference\tinterfering\tbbinteraction" << endl;
    
    for (pair<Residue*,Residue*> cont : _contact_residues) {
        Residue* R_pep = cont.first;
        Residue* R_prot = cont.second;
        mstreal contact = _all_types_contacts["contact"].degree(R_pep,R_prot);
        mstreal interference = _all_types_contacts["interference"].degree(R_pep,R_prot);
        mstreal interfering = _all_types_contacts["interfering"].degree(R_pep,R_prot);
        mstreal bbInteraction = _all_types_contacts["bbinteraction"].degree(R_pep,R_prot);
        out << R_pep->getNum() << "\t" << R_pep->getChainID() << "\t";
        out << R_prot->getNum() << "\t" << R_prot->getChainID() << "\t";
        out << R_pep->getChainID() << R_pep->getNum() << "-" << R_prot->getChainID() << R_prot->getNum() << "\t";
        out << contact << "\t";
        out << interference << "\t";
        out << interfering << "\t";
        out << bbInteraction << endl;
    }
    out.close();
}


void interfaceCoverage::writeAllAlignedSeedsInfo(string outDir) {
    cout << "write all aligned seeds to file..." << endl;
    fstream output;
    string output_path = outDir + "aligned_seeds.tsv";
    MstUtils::openFile(output, output_path, fstream::out);
    //header
    output << "seed_name\tchain_id\tseed_n_terminal_res\tpeptide_n_terminal_res\tlength\tmatch_number\tmatch_rmsd\tsequence_match\trmsd\tavg_angle" << endl;
    //find the segments that contain the peptide residue
    for (int segment_length = 1; segment_length <= chainSubsegmentsBins.size(); segment_length++){
        for (int peptide_position = 0; peptide_position < chainSubsegmentsBins[segment_length-1].size(); peptide_position++) {
            sortedBins& peptide_segment_bins = chainSubsegmentsBins[segment_length-1][peptide_position];
            vector<seedSubstructureInfo> all_seeds = peptide_segment_bins.getAllSeeds();
            
            for (seedSubstructureInfo seed : all_seeds) {
                output << seed.structure_name << "\t";
                output << seed.chain_ID << "\t";
                output << seed.res_idx << "\t";
                output << peptide_position << "\t";
                output << seed.res_length << "\t";
                output << seed.match_number << "\t";
                output << seed.match_rmsd << "\t";
                output << seed.sequence_match << "\t";
                output << seed.rmsd << "\t";
                output << seed.alignment_cos_angle;
                output << endl;
            }
        }
    }
    output.close();
}

void interfaceCoverage::writeBestAlignedSeeds(string outDir, int numSeeds, bool write_structures) {
    cout << "write " << numSeeds << " seeds with lowest RMSD to the peptide to file..." << endl;
    fstream output;
    string output_path = outDir + "best_aligned_seeds.tsv";
    MstUtils::openFile(output, output_path, fstream::out);
    string seedOutDir = "";
    if (write_structures) {
        seedOutDir = outDir + "best_aligned_seeds/";
        MstSys::cmkdir(seedOutDir);
    }
    
    //header
    output << "seed_name\tchain_id\tseed_n_terminal_res\tpeptide_n_terminal_res\tlength\tmatch_number\tmatch_rmsd\tsequence_match\trmsd\tavg_angle\tcontacts" << endl;
    //find the segments that contain the peptide residue
    for (int segment_length = 1; segment_length <= chainSubsegmentsBins.size(); segment_length++){
        for (int peptide_position = 0; peptide_position < chainSubsegmentsBins[segment_length-1].size(); peptide_position++) {
            sortedBins& peptide_segment_bins = chainSubsegmentsBins[segment_length-1][peptide_position];
            vector<seedSubstructureInfo> best_seeds = peptide_segment_bins.getLowestValuePopulatedBin();
            sort(best_seeds.begin(),best_seeds.end());
            
            int avail_seeds = min(int(best_seeds.size()),numSeeds);
            for (int i = 0; i < avail_seeds; i++) {
                seedSubstructureInfo& seed_info = best_seeds[i];
                
                //get the seed segment and write it to a file
                Structure* seed_segment = getSeedSegment(seed_info);
                if (write_structures) seed_segment->writePDB(seedOutDir+MstUtils::toString(segment_length)+"_"+MstUtils::toString(peptide_position)+"_"+seed_segment->getName()+".pdb");
                
                //get the contacts
                set<pair<int,int>> seed_protein_contacts = getContacts(seed_segment->getAtoms(), peptide_position);
                
                //write out contacts in an interpretable format
                stringstream ss;
                int count = 0;
                for (auto cont : seed_protein_contacts) {
                    if (count != 0) ss << " ";
                    ss << cont.first << "," << cont.second;
                    count++;
                }
                
                output << seed_info.structure_name << "\t";
                output << seed_info.chain_ID << "\t";
                output << seed_info.res_idx << "\t";
                output << peptide_position << "\t";
                output << seed_info.res_length << "\t";
                output << seed_info.match_number << "\t";
                output << seed_info.match_rmsd << "\t";
                output << seed_info.sequence_match << "\t";
                output << seed_info.rmsd << "\t";
                output << seed_info.alignment_cos_angle << "\t";
                output << ss.str();
                output << endl;
                
                delete seed_segment;
            }
        }
    }
    output.close();
}

void interfaceCoverage::writeSegmentCoverage(string outDir) {
    //peptide segments (residues)
    fstream out;
    cout << "residue coverage..." << endl;
    string output_path = outDir + "covered_subsegments.tsv";
    MstUtils::openFile(out, output_path, fstream::out);
    
    //write the header
    out << "segment_length\tpeptide_position\tmin_rmsd\tmax_rmsd\tnumber_aligned_seeds" << endl;
    
    for (int segment_length = 1; segment_length <= chainSubsegmentsBins.size(); segment_length++){
        for (int peptide_position = 0; peptide_position < chainSubsegmentsBins[segment_length-1].size(); peptide_position++) {
            vector<tuple<mstreal,mstreal,long>> subsegment_bin_aligned_seeds = chainSubsegmentsBins[segment_length-1][peptide_position].getSeedsByBin();
            for (int bin = 0; bin < subsegment_bin_aligned_seeds.size(); bin++) {
                out << segment_length << "\t";
                out << peptide_position << "\t";
                out << get<0>(subsegment_bin_aligned_seeds[bin]) << "\t";
                out << get<1>(subsegment_bin_aligned_seeds[bin]) << "\t";
                out << get<2>(subsegment_bin_aligned_seeds[bin]) << endl;
            }
        }
    }
    out.close();
}

set<Residue*> interfaceCoverage::getBindingSiteResByDistance(mstreal distanceCutoff) {
    set<Residue*> bindingSiteRes;
    
    // This requires sidechain atoms, so we need to use original complex
    for (Atom* pepA : original_complex->getChainByID(peptide_chain->getID())->getAtoms()) {
        for (Chain* C : target_chains) {
            for (Atom* protA : original_complex->getChainByID(C->getID())->getAtoms()) {
                mstreal distance = pepA->distance(protA);
                if (distance <= distanceCutoff) bindingSiteRes.insert(protA->getResidue());
            }
        }
    }
    return bindingSiteRes;
}

void interfaceCoverage::resetBins() {
    for (int i = 0; i < chainSubsegmentsBins.size(); i++) {
        for (int j = 0; j < chainSubsegmentsBins[i].size(); j++) {
            chainSubsegmentsBins[i][j].reset();
        }
    }
}


// --- Protected --- //

void interfaceCoverage::setParams(string _RL_path) {
    // Import Rotamer Library
    RL_path = _RL_path;
    RL.readRotamerLibrary(RL_path);
    
    seed_chain_id = "0";
    
    //contacts
    cd_threshold = .05;
    int_threshold = .01;
    bbInt_cutoff = 3.25;
    
    //segments
    max_seed_length = 20;
    
    //coverage
    max_rmsd = 2.0; //rmsd to the peptide
    match_rmsd_cutoff = 0.0;  //rmsd of the match when generating seeds (not applied when 0)
    seq_const = false;
    match_number_cutoff = 100000;
}

map<string,contactList> interfaceCoverage::defineContacts(Structure& complex, vector<Residue*> peptide_residues) {
    
    // Structural Element: Residues
    cout << "Total number of residue elements: " << peptide_residues.size() << endl;
    
    // Structural Element: Contacts
    ConFind C(&RL,complex, true);
    
    map<string,contactList> all_types_contacts_map;
    string contact_type;
    contact_type = "contact";
    all_types_contacts_map[contact_type] = generalUtilities::getContactsWith(peptide_residues, C, cd_threshold, contact_type);
    contact_type = "interfering";
    all_types_contacts_map[contact_type] = generalUtilities::getContactsWith(peptide_residues, C, int_threshold, contact_type);
    contact_type = "interfered";
    all_types_contacts_map[contact_type] = generalUtilities::getContactsWith(peptide_residues, C, int_threshold, contact_type);
    contact_type = "bbinteraction";
    all_types_contacts_map[contact_type] = generalUtilities::getContactsWith(peptide_residues, C, bbInt_cutoff, contact_type);
    
    return all_types_contacts_map;
}

void interfaceCoverage::prepareForTERMExtension() {
    vector<Residue*> protein_res;
    map<Residue*,int> complex2target;
    int res_idx = 0;
    for (int chain = 0; chain < complex.chainSize(); chain++) {
        Chain* C = &complex[chain];
        if (C->getID() == peptide_chain->getID()) continue; //get protein only chains
        cout << "Adding chain " << C->getID() << " to new structure" << endl;
        target_chains.push_back(C);
        vector<Residue*> chain_res = C->getResidues();
        protein_res.insert(protein_res.end(),chain_res.begin(),chain_res.end());
        for (Residue* R: chain_res) {
            complex2target[R] = res_idx;
            res_idx++;
        }
    }
    target = new Structure(protein_res);
    target->setName(complex.getName());
    
    cout << "Protein with " << target->chainSize() << " chains and " << target->residueSize() << " residues..." << endl;
    
    vector<Residue*> binding_site;
    for (pair<Residue*,Residue*> cont : contact_residues) {
        if (find(binding_site.begin(),binding_site.end(),cont.second) == binding_site.end()) binding_site.push_back(cont.second);
    }
    
    for (Residue* R : binding_site) bindingSiteRes.push_back(&target->getResidue(complex2target[R]));
}


set<pair<int,int>> interfaceCoverage::getContacts(vector<Atom*> seed_segment, int peptide_position) {
    //construct a new confind object including the protein and the seed
    Structure seed_and_target(*target);
    seed_and_target.addAtoms(seed_segment);
    
    ConFind CD(&RL,seed_and_target);
    
    set<pair<int,int>> contacts;
    //find all contacts between the seed and the protein
    Chain* seed_c = seed_and_target.getChainByID(seed_chain_id);
    vector<Residue*> seed_residues = seed_c->getResidues();
    
    for (Residue* R : seed_residues) {
        contactList R_conts = CD.getContacts(R,cd_threshold);
        for (int i = 0; i < R_conts.size(); i++) {
            /*
             We need three things to get the index of the corresponding peptide residue (in complex)
             The index to the first residue of the peptide (peptide_first_residue_index_in_structure)
             The index (in chain) of the first residue that the seed segment is aligned to (peptide_position)
             The index (in segment) of the current residue (R_conts.residueA(i)->getResidueIndexInChain())
             
             The last assumption is valid because the all atoms in the segment are passed to this function, never a
             subsegment
             */
            Residue* corresponding_peptide_res = &complex.getResidue(peptide_first_residue_index_in_structure + peptide_position + R_conts.residueA(i)->getResidueIndexInChain());
            int seed_r_id = corresponding_peptide_res->getNum();
            int prot_r_id = R_conts.residueB(i)->getNum();
            contacts.emplace(seed_r_id,prot_r_id);
        }
    }
    return contacts;
}

Structure* interfaceCoverage::getSeedSegment(seedSubstructureInfo info) {
    Structure* ext_frag = seeds->getStructureNamed(info.structure_name);
    Chain* C = ext_frag->getChainByID(seed_chain_id);
    vector<Residue*> seed_residues = C->getResidues();
    vector<Residue*> seed_segment_res;
    for (Residue* R : seed_residues) {
        if ((R->getResidueIndexInChain() >= info.res_idx) & (R->getResidueIndexInChain() < info.res_idx + info.res_length)) seed_segment_res.push_back(R);
    }
    Structure* seed_segment = new Structure(seed_segment_res);
    string segment_name = MstUtils::toString(info.rmsd) + "_" + info.structure_name;
    seed_segment->setName(segment_name);
    delete ext_frag;
    return seed_segment;
}

vector<string> pathFromCoveringSeeds::getCoveringPath(int maxSeedLength, int minSeedLength) {
    
    vector<seedSubstructureInfo> pathSeedInfo;
    seedNames.clear();
    coveringSeeds.clear();
    coveringResiduesPath.clear();
    coveringPathString.clear();
    
    /**
     The goal of this function is to find the seeds that fuse to give a backbone that is as representative of the native peptide backbone as
     possible. As our method works by finding overlaps between seeds (which must be at least two residues long) we require that the set
     of covering seeds also contains these overlaps. To achieve this, we cover *segments* of the peptide. For example, by covering all
     three residue segments of the peptide, we guarantee that there is no junction between seeds without at least a two residue overlap.
     
     To obtain the minimum number of low RMSD seeds required to cover the peptide, we employ a greedy algorithm. This algorithm works
     by selecting the seed that covers the most not-yet-covered segments and breaking any ties by the seed with lower RMSD.
     */
    
    set<residueSegment> peptideResidueSegmentsCopy = peptideResidueSegments;
    
    maxSeedLength = min(maxSeedLength,coverage->max_allowable_segment_length);
    
    MstUtils::assert(minSeedLength <= maxSeedLength,"minSeedLength must not be greater than maxSeedLength","pathFromCoveringSeeds::getBestCoveringSeed");
    /*
     The minSeedLength must be at least segmentLength. Smaller seed segments will not be able to
     cover the peptide segments
     */
    minSeedLength = max(minSeedLength,segmentLength);
    
    if (verbose) cout << "Finding covering seeds with minLength: " << minSeedLength << " and maxSeedLength: " << maxSeedLength << endl;
    
    // Continue until all residues have been covered
    int cycle = 0;
    while (!peptideResidueSegments.empty()) {
        cycle++;
        if (verbose) {
            cout << "Cycle " << cycle << endl;
            cout << "peptideResidueSegments: " << peptideResidueSegments.size() << endl;
        }
        
        // get the seed segment that covers the most peptide residues / has lowest RMSD
        pair<seedSubstructureInfo,residueSegment> bestCoveringSeed = getBestCoveringSeed(maxSeedLength,minSeedLength);
        if (bestCoveringSeed.second.empty()) break; //no seeds remaining to cover the peptide, take what is covering and make segments
        if (verbose) cout << "Selected " << bestCoveringSeed.first.structure_name << "_" << bestCoveringSeed.first.res_idx << "_" << bestCoveringSeed.first.res_length << " as the best covering seed aligned to position" << bestCoveringSeed.second.front()->getResidueIndexInChain() << " on the peptide" << endl;
            
        // keep track of seeds that were added, so that only segments from new seeds can be added
        if (forceChimera) seedNames.insert(bestCoveringSeed.first.structure_name);
        auto ret = coveringSeeds.insert(pair<int,seedSubstructureInfo>(bestCoveringSeed.second.front()->getResidueIndexInChain(),bestCoveringSeed.first));
        if (!ret.second) MstUtils::error("Could not insert seed to coveringSeeds, likely a duplicate","pathFromCoveringSeeds::getCoveringPath");
        
        // remove the residues from the set
        bool removeSeg = true;
        set<residueSegment> residueSegmentsCoveredBySeed = getResidueSegmentsCoveredBySeed(bestCoveringSeed.second,removeSeg);
        
        if (verbose) {
            cout << "Covered " << residueSegmentsCoveredBySeed.size() << " peptide residue segments... " << peptideResidueSegments.size() << " segments remaining"  << endl;
            cout << "The following segments were considered covered and removed: " << endl;
            for (residueSegment segment : residueSegmentsCoveredBySeed) {
                cout << *segment.front() << "-" << *segment.back() << endl;
            }
        }
    }
    
    if (coveringSeeds.empty()) {
        cout << "Couldn't find any covering seeds, terminating..." << endl;
        return {};
    }
    cout << "Done finding covering seeds, now try to select residues for fusion" << endl;


    // Get the path residues (as separate vectors if necessary)
    coveringResiduesPath = getPathResidues();
    
    // Reset the coverage elements
    peptideResidueSegments = peptideResidueSegmentsCopy;
    
    // Convert to string (multiple if necessary)
    for (auto coveringRes : coveringResiduesPath) coveringPathString.push_back(PathSampler::getPathFromResidues(coveringRes.second));
    return coveringPathString;
}

void pathFromCoveringSeeds::writeCoveringSeeds(string outputPath) {
    cout << "Now writing seeds covering the peptide: " << endl;
    int count = 0;
    for (auto seed : coveringSeeds) {
        cout << "Seed number: " << count << endl;
        cout << "Seed name: " << seed.second.structure_name << endl;
        cout << "Alignment to peptide:\t" << seed.first << endl;
        cout << "Number of residues:\t" << seed.second.res_length << endl;
        cout << "RMSD: " << seed.second.rmsd << endl;
        
        string seedName = seed.second.structure_name + "-" + MstUtils::toString(seed.first) + "-" + MstUtils::toString(seed.second.res_length) + "-" + MstUtils::toString(seed.second.rmsd);
      
        Structure seedSubstructure = seed.second.loadSeedSubstructureFromCache(seedCache,true);
      
        cout << "Writing pdb for seed substructure with name: " << seedName << endl;
        seedSubstructure.writePDB(outputPath+seedName+".pdb");
        
        count++;
    }
}

pair<seedSubstructureInfo,vector<Residue*>> pathFromCoveringSeeds::getBestCoveringSeed(int maxSeedLength, int minSeedLength) {
    /*
     Count how many not-yet-covered peptide residues would be covered by each seed segment, and take
     the one that would cover the most peptide residues. If there is a tie between seed segments,
     take the one with the lower RMSD. Note that this introduces a slight bias for shorter seeds,
     but I'm assuming that it's not an issue (and this can be controlled with minSeedLength)
     */
    
    multimap<int,pair<sortedBins&,vector<Residue*>>> sortedSeeds;
    for (int i = 0; i < coverage->chainResidueSubsegments.size(); i++) {
        // i+1 is the length of the seeds at this position of chainResidueSubsegments
        if ((i+1 >= minSeedLength) && (i+1 <= maxSeedLength)) {
            //seed length is acceptable
            for (int j = 0; j < coverage->chainResidueSubsegments[i].size(); j++) {
                vector<Residue*> peptideResidues = coverage->chainResidueSubsegments[i][j];
                sortedBins& bins = coverage->chainSubsegmentsBins[i][j];
                if (bins.allBinsEmpty()) continue; //no covering seeds, skip this subsegment
                set<residueSegment> coveredSegments = getResidueSegmentsCoveredBySeed(peptideResidues,false);
                if (coveredSegments.size() > 0) sortedSeeds.emplace(coveredSegments.size(),pair<sortedBins&,vector<Residue*>>(bins,peptideResidues));
            }
        }
    }
    
    if (sortedSeeds.empty()) return pair<seedSubstructureInfo,vector<Residue*>>(seedSubstructureInfo(),{});
    
    // get the highest number of covered residues in the multimap
    int mostCoveredRes = sortedSeeds.rbegin()->first;
    auto ret = sortedSeeds.equal_range(mostCoveredRes);
    
    // get lowest RMSD seed segment
    mstreal lowestRMSD = DBL_MAX;
    pair<seedSubstructureInfo,vector<Residue*>> bestSeed;
    for (auto it = ret.first; it != ret.second; it++) {
        if (!it->second.first.areBinsSorted()) it->second.first.sortBins();
        seedSubstructureInfo seed = it->second.first.getLowestRMSDSeed(seedNames);
        if (seed.structure_name == "") continue;
        
        mstreal rmsd = seed.rmsd;
        if (rmsd < lowestRMSD) {
            lowestRMSD = rmsd;
            bestSeed = pair<seedSubstructureInfo,vector<Residue*>>(seed,it->second.second);
        }
    }
    return bestSeed;
}

set<vector<Residue*>> pathFromCoveringSeeds::getResidueSegmentsCoveredBySeed(vector<Residue*> peptideResiduesCoveredBySeedSegment, bool removeSeg) {
    set<residueSegment> coveredSegments;
    for (int i = 0; i < peptideResiduesCoveredBySeedSegment.size() - segmentLength + 1; i++) {
        residueSegment seedSegment = vector<Residue*>(peptideResiduesCoveredBySeedSegment.begin() + i,peptideResiduesCoveredBySeedSegment.begin() + i + segmentLength);
        
        if (peptideResidueSegments.count(seedSegment) == 1) {
            coveredSegments.insert(seedSegment);
            if (removeSeg) peptideResidueSegments.erase(seedSegment);
        }
    }
    return coveredSegments;
}

map<int,vector<Residue*>> pathFromCoveringSeeds::getPathResidues() {
    /*
     Walk through coveringSeeds, which is sorted by N-terminal residues of the seed segments. Jumping
     points between the seeds are at the center of their overlaps.
     */
    if (coveringSeeds.empty()) {
        cout << "No covering seeds, finding set..." << endl;
        getCoveringPath();
    }
    
    map<int,vector<Residue*>> pathRes;
    vector<Residue*> pathResSegment;
    int currPositionInPeptide = 0, segmentStartPositionInPeptide = 0;
    for (auto it = coveringSeeds.begin(); it != coveringSeeds.end(); it++) {
        /*
         Choose where to start in the seed (positionInSeed) based on where we are in the peptide
         (positionInPeptide) and how the seed maps to the peptide (it->first is the alignment of the
         N-terminal residue of the seed to the peptide)
         
         e.g. if terminalRes = 1
         
         0-1-2-3-4-5-6-7    Peptide Residues
         
         0-1-2-3-4          Seed 1 Residues
                \           crossing point
               0-1-2-3-4    Seed 2 Residues
         
         The starting point within seed 2 is found by taking the position in the peptide (which is 4,
         after adding the first four seed residues) and subtracting the alignment of seed 2 (3).
         
         21/02/03 update: now the class can handle disjoint covering segments. If a region of the
         peptide is not covered (there is no seed aligned where it should be) then the path residue
         vector stops growing, and another one is started. To avoid duplicated residues in the final
         paths, up to the first terminalRes will be truncated from the seed if necessary.
         */
        int positionInSeed = currPositionInPeptide - it->first;
        if (positionInSeed < terminalRes) {
            // we are not where we expected to be in the next seed, there must be a missing seed
            if (it != coveringSeeds.begin()) {
                if (verbose) cout << "adding path segment with " << pathResSegment.size() << " residues" << endl;
                pathRes[segmentStartPositionInPeptide] = pathResSegment;
                pathResSegment.clear();
            }
            // a negative position in seed indicates that we skipped some residues
            // find how many we skipped and set position in seed to beginning
            currPositionInPeptide = currPositionInPeptide - min(positionInSeed,0);
            segmentStartPositionInPeptide = currPositionInPeptide;
            positionInSeed = max(positionInSeed,0);
        }
        
        vector<Residue*> seedChainResidues = it->second.loadSeedSubstructureFromCache(seedCache);
        // if last seed, go to the end
        int skipRes = (next(it) == coveringSeeds.end()) ? 0 : terminalRes;
        for (int i = positionInSeed; i < seedChainResidues.size() - skipRes; i++) {
            pathResSegment.push_back(seedChainResidues[i]);
            currPositionInPeptide++;
        }
    }
    //final push to vector
    if (verbose) cout << "adding path segment with " << pathResSegment.size() << " residues" << endl;
    pathRes[segmentStartPositionInPeptide] = pathResSegment;
        
    return pathRes;
}


mstreal coverageBenchmarkUtils::writeRMSDtoFile(string outputPath, Structure& fusedBackbone, Structure& nativePeptide) {
    //get all of the fused backbone residues and construct a map for searching
    map<int,Residue*> nativePeptideResMap;
    for (Residue* R: nativePeptide.getResidues()) {
        nativePeptideResMap[R->getNum()] = R;
    }
    
    fstream out;
    string outputFile = outputPath+"rmsdByRes.tsv";
    MstUtils::openFile(out,outputFile,fstream::out);
    
    //header line (structure A is fused, structure B is native)
    out << "structureA_residue_id\tstructureA_residue_num\tstructureA_residue_chain_id\tstructureA_residue_aa\tstructureB_residue_id\tstructureB_residue_num\tstructureB_residue_chain_id\tstructureB_residue_aa\trmsd\tcos_angle" << endl;
    
    vector<Atom*> S1_paired_atoms;
    vector<Atom*> S2_paired_atoms;
    
    int pairedPos = 0;
    for (int idx = 0; idx < fusedBackbone.residueSize(); idx++) {
        Residue* R1 = &fusedBackbone.getResidue(idx);
        // see if we can find residue in map
        Residue* R2;
        if (nativePeptideResMap.find(R1->getNum()) != nativePeptideResMap.end()) R2 = nativePeptideResMap[R1->getNum()];
        else continue;
        
        vector<Atom*> R1_atoms = R1->getAtoms();
        vector<Atom*> R2_atoms = R2->getAtoms();
        
        mstreal rmsd = RMSDCalculator::rmsd(R1_atoms,R2_atoms);
        
        // compute the cosine angle between residues
        mstreal cos_angle = generalUtilities::cosAngleBetweenNormalVectors(R1, R2);
        
        out << R1->getResidueIndexInChain() << "\t" << R1->getNum() << "\t" << R1->getChainID() << "\t" << R1->getName() << "\t";
        out << R2->getResidueIndexInChain() << "\t" << R2->getNum() << "\t" << R2->getChainID() << "\t" << R2->getName() << "\t";
        out << rmsd << "\t" << cos_angle << endl;
        
        S1_paired_atoms.insert(S1_paired_atoms.end(),R1_atoms.begin(),R1_atoms.end());
        S2_paired_atoms.insert(S2_paired_atoms.end(),R2_atoms.begin(),R2_atoms.end());
        
        pairedPos++;
    }
    cout << "found " << pairedPos << " paired positions between provided structures" << endl;
    return RMSDCalculator::rmsd(S1_paired_atoms,S2_paired_atoms);
}

void coverageBenchmarkUtils::writeContactstoFile(string outputPath, interfaceCoverage *IC, Structure& fusedPathandTarget, set<string> peptideChains, string rotLibFile) {
    vector<Residue*> fusedPathResidues;
    for (Residue *R : fusedPathandTarget.getResidues()) if (peptideChains.find(R->getChainID()) != peptideChains.end()) fusedPathResidues.push_back(R);
    
    // Get contacts between fused path residues and target
    map<string, contactList> all_types_contacts;
    all_types_contacts = IC->defineContacts(fusedPathandTarget,fusedPathResidues);
    
    set<pair<Residue*,Residue*>> contact_residues;
    contact_residues = generalUtilities::mergeContactLists({all_types_contacts["contact"],all_types_contacts["interfering"],all_types_contacts["interfered"],all_types_contacts["bbinteraction"]});
    
    cout << "Total number of contacts between fused path and protein " << contact_residues.size() << endl;
    
    IC->writeContacts(outputPath, contact_residues, all_types_contacts);
}

void coverageBenchmarkUtils::getResiduesFromMap(string resMapPath, Structure& structureA, Structure& structureB, vector<Residue*>& selectedResA, vector<Residue*>& selectedResB) {
    map<string,Chain*> chainMapA, chainMapB;
    for (int i = 0; i < structureA.chainSize(); i++) chainMapA[structureA.getChain(i).getID()] = &structureA.getChain(i);
    for (int i = 0; i < structureB.chainSize(); i++) chainMapB[structureB.getChain(i).getID()] = &structureB.getChain(i);
    map<pair<Chain*,int>,Residue*> resMapA, resMapB;
    for (Residue* R : structureA.getResidues()) resMapA[pair<Chain*,int>(R->getChain(),R->getNum())] = R;
    for (Residue* R : structureB.getResidues()) resMapB[pair<Chain*,int>(R->getChain(),R->getNum())] = R;
    
    // open the residue map
    vector<string> lines = MstUtils::fileToArray(resMapPath);
    bool header = true;
    for (string line : lines) {
        if (header == true) {
            header = false;
            continue;
        }
        vector<string> lineSplit = MstUtils::split(line,"\t");
        
        // this file is expected to be generated by writeRMSDtoFile()
        string chainA = lineSplit[2];
        int resNumA = MstUtils::toInt(lineSplit[1]);
        
        //cout << chainA << " " << resNumA << " " << endl;
        
        Residue* resA = resMapA.at(pair<Chain*,int>(chainMapA.at(chainA),resNumA));
        selectedResA.push_back(resA);
        
        string chainB = lineSplit[6];
        int resNumB = MstUtils::toInt(lineSplit[5]);
        
        //cout << chainB << " " << resNumB << " " << endl;
        
        Residue* resB = resMapB.at(pair<Chain*,int>(chainMapB.at(chainB),resNumB));
        selectedResB.push_back(resB);
    }
}

void coverageBenchmarkUtils::fuseCoveringSeeds(interfaceCoverage* IC, bool force_chimera, int max_seed_length_fuse, string fusDir, string pdb_id, Structure& target, Structure& complex, bool two_step_fuse, string RL) {
    int segmentLength = 3;
    pathFromCoveringSeeds generatePath(IC,segmentLength,force_chimera);

    cout << "Find covering residues and generate a path" << endl;
    vector<string> pathString = generatePath.getCoveringPath(max_seed_length_fuse);
    if (pathString.empty()) {
        cout << "There were no paths found" << endl;
        return;
    }
    else {
        cout << pathString.size() << " path(s) were generated when trying to cover the peptide: " << endl;
        for (string path : pathString) {
            cout << path << endl;
        }
    }

    generatePath.writeCoveringSeeds(fusDir+"_");
    
    /**
     Now fuse each set of seeds together and write the output to a pdb file
     */
    Structure fusedPath, fusedPathAndTarget(target);
    Structure peptideStructure(*IC->getPeptideChain());
    
    int overlapLength = 1;
    PathSampler sampler(&target,overlapLength); int count = 0;
    sampler.setTwoStepFuse(two_step_fuse);
    set<string> fusedPathChains; string fusedPeptideChainsString;
    for (auto it : generatePath.getCoveringResiduePath()) {
        int peptideStartPosition = it.first;
        vector<Residue*> pathResidues = it.second;
        cout << "Fuse path: " << count << endl;
        vector<PathResult> results;
        bool ignore_clashes = true;
        bool clashes = sampler.emplacePathFromResidues(pathResidues, results, {}, ignore_clashes);
        cout << "Does the fused output clash?: " << clashes << endl;
        PathResult fusedCoveringPath = results.back();

        fusionOutput fuseOutput = fusedCoveringPath.getFuserScore();
        cout << "Printing fuser output... " << endl;
        cout << "Bond score: " << fuseOutput.getBondScore() << endl;
        cout << "Angle score: " << fuseOutput.getAngleScore() << endl;
        cout << "Dihedral score: " << fuseOutput.getDihedralScore() << endl;
        cout << "RMSD score: " << fuseOutput.getRMSDScore() << endl;
        cout << "TotRMSDScore: " << fuseOutput.getTotRMSDScore() << endl;
        cout << "Score: " << fuseOutput.getScore() << endl;
        
        // get fused path and renumber the residues according to the native peptide
        Structure fusedPathSection;
        fusedCoveringPath.getFusedPathOnly(fusedPathSection);
        
        cout << "Fused path with " << fusedPathSection.residueSize() << " residues" << endl;
        
        // set number to match the peptide and set residue name to "UNK"
        int segmentPos = 0;
        for (Residue* R : fusedPathSection.getResidues()) {
            R->setNum(complex.getChainByID(IC->getPeptideChain()->getID())->getResidue(peptideStartPosition+segmentPos).getNum());
            R->setName("UNK");
            segmentPos++;
        }

        Chain* fusedPathChain = new Chain(fusedPathSection.getChain(0));
        fusedPathChain->setID(MstUtils::toString(count));
        fusedPath.appendChain(fusedPathChain);
        Chain* fusedPathChainCopy = new Chain(*fusedPathChain);
        fusedPathAndTarget.appendChain(fusedPathChainCopy);
        
        fusedPathChains.insert(fusedPathChain->getID());
        fusedPeptideChainsString += fusedPathChain->getID();
        
        count++;
    }
    
    // add peptide/protein chains to structure name
    string proteinChainsString;
    for (int i = 0; i < target.chainSize(); i++) proteinChainsString += target.getChain(i).getID();
    string pdb_name = pdb_id + "_" + proteinChainsString + "_" + fusedPeptideChainsString;
    
    mstreal rmsd = coverageBenchmarkUtils::writeRMSDtoFile(fusDir+"_"+pdb_name+"_",fusedPath,peptideStructure);
    cout << "The RMSD between the peptide and the fused backbone is: " << rmsd << endl;
    
    coverageBenchmarkUtils::writeContactstoFile(fusDir+"_"+pdb_name+"_", IC, fusedPathAndTarget, fusedPathChains, RL);
    
    fusedPath.writePDB(fusDir+"_"+pdb_name+"_pathOnly.pdb");
    fusedPathAndTarget.writePDB(fusDir+"_"+pdb_name+".pdb");
}

