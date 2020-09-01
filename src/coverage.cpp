//
//  compareterms.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 1/29/19.
//

#include "coverage.h"

using namespace MST;

/* --------- sortedBins --------- */
void sortedBins::insert(seedSubstructureInfo& entry, mstreal val) {
    if (val > max_val) return;
    if (bins.empty()) MstUtils::error("Tried to insert an entry to an improperly constructed instance of sortedBins");
    bins[val2Bin(val)].push_back(entry);
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
    complex = Structure(RotamerLibrary::getBackbone(S));
    complex.setName(S.getName());
    
    peptide_chain = complex.getChainByID(p_cid);
    peptide_first_residue_index_in_structure = peptide_chain->getResidue(0).getResidueIndex();
    
    //define coverage elements
    cout << "Identifying structural elements..." << endl;
    defineCoverageElements();
    
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
}

void interfaceCoverage::setSeeds(string _binFilePath) {
    if (binFilePath != "") {
        resetBins();
        delete seeds;
    }
    binFilePath = _binFilePath;
    seeds = new StructuresBinaryFile(binFilePath);
    seeds->scanFilePositions();
    seeds->reset();
}

void interfaceCoverage::findCoveringSeeds() {
    resetBins();
    while (seeds->hasNext()) {
        Structure* extended_fragment = seeds->next();
        
        //check if structure meets criteria
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
        mapSeedToChainSubsegments(getBackboneAtoms(seed_C),seed_C->getResidues(),match_rmsd,sequence_match);
        delete extended_fragment;
    }
    seeds->reset();
}

void interfaceCoverage::mapSeedToChainSubsegments(vector<Atom*> seed_atoms, vector<Residue*> seed_residues, mstreal match_rmsd, bool sequence_match) {
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
                if (rmsd < max_rmsd) {
                    Atom* A = seed_segment[0];
                    string structure_name = A->getStructure()->getName();
                    string chain_ID = A->getChain()->getID();
                    
                    //compute the mean cosine angle between residue normal vectors
                    vector<Residue*> seed_residue_segment(seed_residues.begin()+seed_position,seed_residues.begin()+seed_position+segment_length);
                    mstreal alignment_cos_angle = generalUtilities::avgCosAngleBetweenSegments(seed_residue_segment,chainResidueSubsegments[segment_length-1][peptide_position]);
                    
                    seedSubstructureInfo info(structure_name,chain_ID,seed_position,segment_length,match_rmsd,sequence_match,rmsd,alignment_cos_angle);
                    
                    chainSubsegmentsBins[segment_length-1][peptide_position].insert(info, rmsd);
                }
            }
        }
    }
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

void interfaceCoverage::writeContacts(string outDir) {
    //contacts
    cout << "contact info..." << endl;
    fstream out;
    string output_path = outDir + "contacts.tsv";
    MstUtils::openFile(out, output_path, fstream::out);
    
    out << "peptide_residue_number\tpeptide_chain\tprotein_residue_number\tprotein_chain\tcontact_name\tcontact\tinterference\tinterfering\tbbinteraction" << endl;
    
    for (pair<Residue*,Residue*> cont : contact_residues) {
        Residue* R_pep = cont.first;
        Residue* R_prot = cont.second;
        mstreal contact = all_types_contacts["contact"].degree(R_pep,R_prot);
        mstreal interference = all_types_contacts["interference"].degree(R_pep,R_prot);
        mstreal interfering = all_types_contacts["interfering"].degree(R_pep,R_prot);
        mstreal bbInteraction = all_types_contacts["bbinteraction"].degree(R_pep,R_prot);
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
    output << "seed_name\tchain_id\tseed_n_terminal_res\tpeptide_n_terminal_res\tlength\tmatch_rmsd\tsequence_match\trmsd\tavg_angle" << endl;
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

void interfaceCoverage::writeBestAlignedSeeds(string outDir, int numSeeds) {
    cout << "write " << numSeeds << " seeds with lowest RMSD to the peptide to file..." << endl;
    fstream output;
    string output_path = outDir + "best_aligned_seeds.tsv";
    MstUtils::openFile(output, output_path, fstream::out);
    string seedOutDir = outDir + "best_aligned_seeds/";
    MstSys::cmkdir(seedOutDir);
    
    //header
    output << "seed_name\tchain_id\tseed_n_terminal_res\tpeptide_n_terminal_res\tlength\tmatch_rmsd\tsequence_match\trmsd\tavg_angle\tcontacts" << endl;
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
//                seed_segment->writePDB(seedOutDir+seed_segment->getName()+".pdb");
                
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
    max_seed_length = 10;
    
    //coverage
    max_rmsd = 3.0; //rmsd to the peptide
    match_rmsd_cutoff = 0.0;  //rmsd of the match when generating seeds (not applied when 0)
    seq_const = false;
}

void interfaceCoverage::defineCoverageElements() {
    // Structural Element: Residues
    peptide_residues = peptide_chain->getResidues();
    cout << "Total number of residue elements: " << peptide_residues.size() << endl;
    
    // Structural Element: Contacts
    ConFind C(&RL,complex, true);
    
    string contact_type;
    contact_type = "contact";
    all_types_contacts[contact_type] = generalUtilities::getContactsWith(peptide_residues, C, cd_threshold, contact_type);
    contact_type = "interfering";
    all_types_contacts[contact_type] = generalUtilities::getContactsWith(peptide_residues, C, int_threshold, contact_type);
    contact_type = "interfered";
    all_types_contacts[contact_type] = generalUtilities::getContactsWith(peptide_residues, C, int_threshold, contact_type);
    contact_type = "bbinteraction";
    all_types_contacts[contact_type] = generalUtilities::getContactsWith(peptide_residues, C, bbInt_cutoff, contact_type);
    
    
    contact_residues = generalUtilities::mergeContactLists({all_types_contacts["contact"],all_types_contacts["interfering"],all_types_contacts["interfered"],all_types_contacts["bbinteraction"]});
    
    cout << "Total number of contact residue elements: " << contact_residues.size() << endl;
    
    for (pair<Residue*,Residue*> cont : contact_residues) cout << cont.first->getChainID() << cont.first->getNum() << "-" << cont.second->getChainID() << cont.second->getNum() << " ";
    cout << endl;
}

void interfaceCoverage::prepareForTERMExtension() {
    vector<Residue*> protein_res;
    map<Residue*,int> complex2target;
    int res_idx = 0;
    for (int chain = 0; chain < complex.chainSize(); chain++) {
        Chain* C = &complex[chain];
        if (C->getID() == peptide_chain->getID()) continue; //get protein only chains
        cout << "Adding chain " << C->getID() << " to new structure" << endl;
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
    string segment_name = info.structure_name + ":" + MstUtils::toString(info.res_idx)+ "-" + MstUtils::toString(info.res_idx+info.res_length);
    seed_segment->setName(segment_name);
    delete ext_frag;
    return seed_segment;
}


///* --------- benchmarkUtils --------- */
//seedScoring::seedScoring(const string& configFile):dTERMen(configFile){
//  //copy over necessary private variables from dTERMen
//  RL_p = this->getRotamerLibrary();
//  globalAlph = this->getGlobalAlph();
//  
//  //make aa triple letter alphabet
//  for (int k = 0; k < globalAlphabetSize(); k++) triple_letter_alph.push_back(SeqTools::idxToTriple(globalAlph[k]));
//}
//
//vector<mstreal> seedScoring::seedDist(Residue* cen_R) {
//  vector<mstreal> seed_dist;
//  
//  Structure* S = cen_R->getStructure();
//  ConFind C(RL_p, *S, true);
//  
//  // all residues that contact variable positions, but are not variable, are fixed
//  set<Residue*> fixed_set;
//  vector<pair<Residue*, Residue*>> conts = getContactsWith({cen_R}, C, 0); //type = 0 because only "inter" contacts are desired
//  for (int i = 0; i < conts.size(); i++) {
//    Residue* res = conts[i].second;
//    fixed_set.insert(res);
//  }
//  
//  // make sure amino acids at fixed positions are legal
//  for (auto it = fixed_set.begin(); it != fixed_set.end(); ++it) {
//    Residue* res = *it;
//    if (!aaIndexKnown(aaToIndex(res->getName()))) MstUtils::error("fixed position " + MstUtils::toString(*(*it)) + " occupied with unknown amino acid", "seedScoring::dTERMen::buildEnergyTable");
//  }
//  
//  //// compute self energies
//  cout << "computing self energy..." << endl;
//  //print res names?
//  vector<mstreal> selfE = selfEnergies(cen_R, C, true);
//  
//  // compute pair energies for fixed-contacts and add to self energies
//  int i = 0;
//  for (auto it = fixed_set.begin(); it != fixed_set.end(); it++, i++) {
//    Residue* resA = cen_R;
//    Residue* resB = *it;
//    cout << "computing pair energy for positions " << *resA << " x " << *resB << ", " << i+1 << "/" << conts.size() << endl;
//    vector<vector<mstreal>> pairE = pairEnergies(resA, resB, true);
//    
//    // interaction between a variable and a fixed position (goes into self)
//    int bi = aaToIndex(resB->getName());
//    int naa = globalAlphabetSize();
//    for (int a = 0; a < naa; a++) {
//      int ai = aaToIndex(triple_letter_alph[a]);
//      selfE[a] = selfE[a] + pairE[ai][bi];
//      //curiousity
//      cout << "are the index from the loop and the index from the function the same?" << endl;
//      cout << "aaToIndex() -> " << MstUtils::toString(ai) << " index from loop: " << MstUtils::toString(a) << endl;
//    }
//  }
//  return seed_dist;
//}
//
//vector<mstreal> seedScoring::backgroundDist(Residue* R) {
//  
//  // Background energy estimated by counting each amino acid in the fasstDB
//  CartesianPoint backE(globalAlphabetSize(), 0.0);
//  for (int aai = 0; aai < backE.size(); aai++) {
//    backE[aai] = backEner(aai);
//  }
//  
//  //convert energies to probabilities
//  vector<mstreal> p(backE);
//  dTERMen::enerToProb(p);
//  
//  return p;
//}
//
//
//vector<mstreal> seedScoring::uniformDist() {
//  
//  int naa = globalAlphabetSize();
//  vector<mstreal> p(naa, 1.00/naa);
//  
//  return p;
//}

