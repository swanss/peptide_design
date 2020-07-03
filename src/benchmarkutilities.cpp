//
//  seedutilities.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 4/30/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include "benchmarkutilities.h"

/* --------- histogram --------- */

void histogram::readHistFile(string hist_file) {
    cout << "reading histogram from file" << endl;
    bins.clear();
    
    fstream info_file;
    MstUtils::openFile(info_file,hist_file,fstream::in);
    
    string line;
    vector<int> bins_count;
    int line_count = 0;
    while (getline(info_file,line)) {
        line_count++;
        if (line_count == 1) {
            if (line != "lower_bound,upper_bound,value") MstUtils::error("Wrong file type (header does not match)","rejectionSampler::rejectionSampler()");
            continue;
        }
        vector<string> split = MstUtils::split(line,",");
        MstUtils::assert(split.size() == 3); //should have three entries
        if (line_count == 2) {
            min_value = MstUtils::toReal(split[0]);
        }
        bins.push_back(MstUtils::toReal(split[2]));
        max_value = MstUtils::toReal(split[1]); //final value will be the max_value
    }
    bin_size = (max_value - min_value) / bins.size();
    if (line_count == 1) MstUtils::error("File had no entries","rejectionSampler::rejectionSampler()");
    info_file.close();
    cout << "min: " << min_value << "\tmax: " << max_value << "\tbin size: " << bin_size << endl;
}

void histogram::writeHistFile(string hist_file) {
    fstream info_file;
    MstUtils::openFile(info_file,hist_file,fstream::out);
    
    info_file << "lower_bound,upper_bound,value" << endl;
    
    mstreal lower_bound = min_value;
    mstreal upper_bound = min_value + bin_size;
    for (mstreal value : bins) {
        info_file << lower_bound << "," << upper_bound << "," << value << endl;
        lower_bound += bin_size;
        upper_bound += bin_size;
    }
    info_file.close();
}

/* --------- seedCentroidDistance --------- */

seedCentroidDistance::seedCentroidDistance(string list, mstreal _min_value, mstreal _max_value, int num_bins) : min_value(_min_value), max_value(_max_value) {
    bin_counts.resize(num_bins);
    fill(bin_counts.begin(),bin_counts.end(),0);
    
    bin_size = (max_value - min_value) / num_bins;
    cout << "bin size: " << bin_size << endl;
    
    hist.min_value = min_value;
    hist.max_value = min_value + (bin_size*num_bins);
    hist.bin_size = bin_size;
    
    //get the list of binary files/structures
    //e.g. path/to/binary path/to/structure
    vector<string> files = MstUtils::fileToArray(list);
    cout << "there are " << files.size() << " sets of seeds" << endl;
    
    if (files.size() * sample_n > INT_MAX) MstUtils::error("Note: it is possible that the counts in a single bin could exceed the maximum possible value of an 'int'. If this occurs, consider sampling less seeds or refactoring");
    
    for (string file_name : files) {
        cout << file_name << endl;
        vector<string> split = MstUtils::split(file_name);
        string bin_path = split[0];
        string struct_path = split[1];
        
        //get the peptide chain name
        split = MstUtils::split(struct_path,"_");
        string p_id = split[split.size() - 3]; //should get the peptide chain id
        
        cout << "pdb path: " << struct_path << " bin path: " << bin_path << " peptide chain ID: " << p_id << endl;
        
        //sample seeds from binary file and record distance from protein
        Structure complex(struct_path);
        seedStatistics stat(complex,p_id);
        
        StructuresBinaryFile bin(bin_path,true);
        size_t num_seeds = bin.structureCount();
        bin.reset();
        
        mstreal sample_prob = min(1.0,mstreal(sample_n)/mstreal(num_seeds)); //such that, on avg, sample_n seeds are sampled
        Structure* extfrag;
        while (bin.hasNext()) {
            if (MstUtils::randUnit() < sample_prob) extfrag = bin.next();
            else {
                bin.skip();
                continue;
            }
            Chain* C = extfrag->getChainByID("0");
            Structure* seed = new Structure(*C);
            mstreal distance = stat.centroid2NearestProteinAtom(seed);
            
            bin_counts[getBin(distance)] += 1;
            
            delete extfrag;
            delete seed;
        }
    }
    
    //build histogram (normalized by bin with highest count)
    int max_value = 0;
    for (int val : bin_counts) if (val > max_value) max_value = val;
    for (int val : bin_counts) hist.bins.push_back(mstreal(val)/mstreal(max_value));
}


/* --------- rejectionSampler --------- */

rejectionSampler::rejectionSampler(string hist_file) {
    hist.readHistFile(hist_file);
    cout << "loaded the following histogram: " << endl;
    hist.printInfo();
}

bool rejectionSampler::accept(mstreal value) {
    mstreal accept_prob = getVal(value);
    mstreal sampled_value = MstUtils::randUnit();
    if (sampled_value <= accept_prob) return true;
    else return false;
}

mstreal rejectionSampler::getVal(mstreal value) {
    if ((value < hist.min_value) || (value > hist.max_value)) {
        string value_str = MstUtils::toString(value);
        string min_str = MstUtils::toString(hist.min_value);
        string max_str = MstUtils::toString(hist.max_value);
        MstUtils::error("Provided value ("+value_str+") not in rejection sampler range: ["+min_str+","+max_str+"]","rejectionSampler::getCount()");
    }
    int bin_id = floor((value - hist.min_value) / hist.bin_size);
    return hist.bins[bin_id];
}

/* --------- generateRandomSeed --------- */

generateRandomSeed::generateRandomSeed(const string& dbFile, int _max_len) : classifier() {
    cout << "Constructing generateRandomSeed class for sampling segments from the DB" << endl;
    readDatabase(dbFile,1);
    max_len = _max_len;
    
    //build windows
    for (int l = 1; l <= max_len; l++) {
        windows.insert(map<int,vector<pair<int,int>>>::value_type(l,vector<pair<int,int>>()));
        int target_num = numTargets();
        for (int target_id = 0; target_id < target_num; target_id++) {
            Structure* target = getTarget(target_id);
            int chain_size = target->chainSize();
            int chain_res_id = 0; //the index of the N-terminal residue of the current chain
            //only build windows within chains (which are contiguous)
            for (int chain = 0; chain < chain_size; chain++) {
                int window_num = target->getChain(chain).residueSize() - l + 1;
                //add all windows to the vector
                for (int window_id = 0; window_id < window_num; window_id++) {
                    windows[l].push_back(make_pair(target_id,chain_res_id+window_id));
                }
                chain_res_id += target->getChain(chain).residueSize();
            }
        }
    }
};

pair<Structure*,string> generateRandomSeed::getSeed(int len) {
    if (windows.find(len) == windows.end()) MstUtils::error("Seed length must be in [1,max_length] ("+MstUtils::toString(max_len)+")");
    //sample a random window
    int window_id = MstUtils::randInt(windows[len].size());
    int target_id = windows[len][window_id].first;
    int nterm_res_id = windows[len][window_id].second;
    
    //copy the residues from the target protein
    Structure* target = getTarget(target_id);
    Structure* seed = new Structure;
    Chain* C = seed->appendChain("0",false);
    vector<int> res_idx; //for getting the secondary structure in the next section
    for (int i = 0; i < len; i++) {
        res_idx.push_back(nterm_res_id+i);
        Residue* R = new Residue(target->getResidue(res_idx[i]));
        C->appendResidue(R);
    }
    
    //get secondary structure
    string sec_structure = classifier.classifyResInStruct(target,res_idx);
    
    //name the new seed
    //targetid_resid_length
    string seed_name = MstUtils::toString(target_id) + "_" + MstUtils::toString(nterm_res_id) + "_" + MstUtils::toString(len);
    seed->setName(seed_name);
    
    return make_pair(seed,sec_structure);
}

/* --------- structureBoundingBox --------- */

structureBoundingBox::structureBoundingBox(Chain * C, mstreal _pad) : pad(_pad){
    Structure Chain_structure = *C;
    AtomPointerVector structure_atoms = RotamerLibrary::getBackbone(Chain_structure);
    construct_structureBoundingBox(structure_atoms);
}

structureBoundingBox::structureBoundingBox(vector<Residue*> residues, mstreal _pad) : pad(_pad){
    vector<Atom*> atoms;
    for (Residue* R: residues) {
        if (!RotamerLibrary::hasFullBackbone(R)) MstUtils::error("Provided residue is missing backbone atoms");
        vector<Atom*> R_bb_atoms = RotamerLibrary::getBackbone(R);
        atoms.insert(atoms.end(),R_bb_atoms.begin(),R_bb_atoms.end());
    }
    construct_structureBoundingBox(atoms);
}

void structureBoundingBox::construct_structureBoundingBox(AtomPointerVector atoms) {
    cout << "Trying to construct bounding box containing all seed centroids..." << endl;
    xlo = INFINITY;
    ylo = INFINITY;
    zlo = INFINITY;
    xhi = -INFINITY;
    yhi = -INFINITY;
    zhi = -INFINITY;
    
    cout << "Backbone with " << atoms.size() << " number atoms" << endl;
    for (Atom* A : atoms) {
        if (A->getX() < xlo) xlo = A->getX();
        if (A->getX() > xhi) xhi = A->getX();
        if (A->getY() < ylo) ylo = A->getY();
        if (A->getY() > yhi) yhi = A->getY();
        if (A->getZ() < zlo) zlo = A->getZ();
        if (A->getZ() > zhi) zhi = A->getZ();
    }
    
    //add the pad to each boundary
    xlo -= pad; xhi += pad;
    ylo -= pad; yhi += pad;
    zlo -= pad; zhi += pad;
    
    cout << "Bounding box complete with the following boundaries: ";
    cout << "Xlo: " << xlo << ", Xhi: " << xhi << "\t";
    cout << "Ylo: " << ylo << ", Yhi: " << yhi << "\t";
    cout << "Zlo: " << zlo << ", Zhi: " << zhi << "\t";
    cout << endl;
}

/* --------- naiveSeedsfromBin --------- */
naiveSeedsFromBin::naiveSeedsFromBin(Structure& S, string p_id, string seedBinaryPath_in, string sampler_path, mstreal _distance, int _neighbors) : complex(S), seeds(seedBinaryPath_in), sampler(sampler_path), stat(S,p_id) {
    //get target structure
    peptide = complex.getChainByID(p_id);
    vector<Residue*> target_residues;
    for (Residue* R : complex.getResidues()) if (R->getChainID() != p_id) target_residues.push_back(R);
    target = Structure(target_residues);
    
    //generate proximity search for target backbone
    cout << "Extracting target backbone to construct a promixity search object..." << endl;
    if (!RotamerLibrary::hasFullBackbone(target)) cout << "warning: target structure is missing backbone atoms!" << endl;
    target_BB_atoms = RotamerLibrary::getBackbone(target);
    target_BB_structure = Structure(target_BB_atoms);
    target_PS = new ProximitySearch(target_BB_atoms, vdwRadii::maxSumRadii()/2);
    
    //generate proximity search for the seed centroids
    cout << "Loading seeds and calculating centroids to construct a promixity search object..." << endl;
    
    //get file positions and properties
    seeds.scanFilePositions();
    seeds.reset();
    int_properties = seeds.getPropertyNamesInt();
    real_properties = seeds.getPropertyNamesReal();
    
    distance = _distance;
    //  seed_PS = new ProximitySearch(seed_centroids,distance/2);
    neighbors = _neighbors;
    max_attempts = 50;
}

void naiveSeedsFromBin::newPose(string output_path, string out_name, bool position, bool orientation, vector<Residue*> binding_site) {
    
    if ((position == false) && (orientation == false)) MstUtils::error("Neither position or orientation have been selected to be randomized. Nothing to do!");
    cout << "Will randomize: ";
    if (position) cout << "position\t";
    if (orientation) cout << "orientation";
    cout << endl;
    
    //construct bounding box
    structureBoundingBox bounding_box = structureBoundingBox(peptide); //only use if position is randomized
    
    //prepare for writing new seed binary file
    string seedBinaryPath_out = output_path + out_name + ".bin";
    
    StructuresBinaryFile seeds_out(seedBinaryPath_out,false);
    
    //open seed data file
    string seedRetry_out = output_path + out_name + "_newpose_attempts.out";
    fstream retry_out;
    MstUtils::openFile(retry_out, seedRetry_out, fstream::out, "naiveSeedsFromBin::newPose");
    //header
    retry_out << "name\tattempts" << endl;
    
    //randomize seeds and write to file
    int count = 0;
    seeds.reset();
    while (seeds.hasNext()) {
        Structure* extended_fragment = seeds.next();
        Chain* seed_C = extended_fragment->getChainByID(seed_chain_id);
        Structure* seed_structure = new Structure(*seed_C);
        CartesianPoint seed_centroid = AtomPointerVector(seed_structure->getAtoms()).getGeometricCenter();
        
        //randomize position/orientation
        int attempts = transform(seed_structure, bounding_box, position, orientation, seed_centroid);
        
        //write transformed seed to binary file (including meta-data)
        //structure
        seeds_out.appendStructure(seed_structure);
        //meta-data
        for (string prop_name : int_properties) {
            int value = seeds.getStructurePropertyInt(prop_name,extended_fragment->getName());
            seeds_out.appendStructurePropertyInt(prop_name,value);
        }
        for (string prop_name : real_properties) {
            mstreal value = seeds.getStructurePropertyReal(prop_name,extended_fragment->getName());
            seeds_out.appendStructurePropertyReal(prop_name,value);
        }
        
        //write number of attempts to info file
        retry_out << seed_structure->getName() << "\t" << attempts << endl;
        
        delete extended_fragment;
        delete seed_structure;
        count++;
    }
    seeds.reset();
    int num_seeds = seeds.structureCount();
    cout << "There are " << num_seeds << " seeds in the input binary file. After randomization, there are " << count << " seeds in the output binary file" << endl;
    retry_out.close();
}


int naiveSeedsFromBin::transform(Structure* seed, structureBoundingBox& bounding_box, bool position, bool orientation, CartesianPoint new_centroid) {
    //make a copy of the atoms in the structure, in case multiple transformations must be sampled
    vector<Atom*> original_seed_atoms = seed->getAtoms(); vector<Atom*> original_seed_atoms_copy;
    for (Atom* A : original_seed_atoms) {
        Atom* A_copy = new Atom(A);
        original_seed_atoms_copy.push_back(A_copy);
    }
    
    TransformFactory tf;
    //initial transform to bring the seed centroid to the origin
    CartesianPoint seed_original_centroid = AtomPointerVector(original_seed_atoms_copy).getGeometricCenter();
    CartesianPoint centroid_translate = CartesianPoint(0,0,0) - seed_original_centroid; //this vector should cancel with seed_original_centroid
    Transform translate_to_origin = tf.translate(centroid_translate);
    
    //construct all other transforms
    //only used if orientation is selected
    Transform rotate_around_x, rotate_around_y, rotate_around_z;
    
    //only used if translation is selected
    Transform translate_to_bounding_box_origin;
    Transform translate_into_box;
    
    //only used if translation is not selected
    Transform translate_to_centroid;
    
    //final transform
    Transform sample_new_pose;
    
    int attempts = 0;
    bool seed_clash = true;
    while (seed_clash) {
        if (attempts > 0) {
            generalUtilities::copyAtomCoordinates(seed, original_seed_atoms_copy);
        }
        attempts++;
        
        if (orientation) {
            //generate a random rotation
            mstreal x_angle = MstUtils::randUnit() * 360;
            mstreal y_angle = MstUtils::randUnit() * 360;
            mstreal z_angle = MstUtils::randUnit() * 360;
            rotate_around_x = tf.rotateAroundX(x_angle);
            rotate_around_y = tf.rotateAroundY(y_angle);
            rotate_around_z = tf.rotateAroundZ(z_angle);
        }
        
        //generate the translation
        if (position) {
            //      // sample and check if near seeds
            //      mstreal x_pos;
            //      mstreal y_pos;
            //      mstreal z_pos;
            //      bool near_seeds = false;
            //      while (!near_seeds) {
            //        x_pos = bounding_box.xlo + (MstUtils::randUnit() * (bounding_box.xhi - bounding_box.xlo));
            //        y_pos = bounding_box.ylo + (MstUtils::randUnit() * (bounding_box.yhi - bounding_box.ylo));
            //        z_pos = bounding_box.zlo + (MstUtils::randUnit() * (bounding_box.zhi - bounding_box.zlo));
            //        CartesianPoint new_centroid(x_pos,y_pos,z_pos);
            //        vector<int> near_points = seed_PS->getPointsWithin(new_centroid, 0.0, distance);
            //        if (near_points.size() > neighbors) near_seeds = true;
            //      }
            // sample until position is accepted
            mstreal x_pos;
            mstreal y_pos;
            mstreal z_pos;
            bool accept = false;
            while (!accept) {
                x_pos = bounding_box.xlo + (MstUtils::randUnit() * (bounding_box.xhi - bounding_box.xlo));
                y_pos = bounding_box.ylo + (MstUtils::randUnit() * (bounding_box.yhi - bounding_box.ylo));
                z_pos = bounding_box.zlo + (MstUtils::randUnit() * (bounding_box.zhi - bounding_box.zlo));
                CartesianPoint new_centroid(x_pos,y_pos,z_pos);
                mstreal distance = stat.point2NearestProteinAtom(new_centroid);
                accept = sampler.accept(distance);
            }
            
            translate_into_box = tf.translate(x_pos,y_pos,z_pos);
        } else {
            //to the provided centroid
            translate_to_centroid = tf.translate(new_centroid);
        }
        
        //combine all transformations
        if (position && orientation) {
            sample_new_pose = translate_into_box * translate_to_bounding_box_origin * rotate_around_z * rotate_around_y * rotate_around_x * translate_to_origin; //matrix multiplication progresses right to left.... right? lol
        } else if (!position && orientation) {
            sample_new_pose = translate_to_centroid * rotate_around_z * rotate_around_y * rotate_around_x * translate_to_origin;
        } else if (position && !orientation) {
            sample_new_pose = translate_into_box * translate_to_bounding_box_origin * translate_to_origin;
        }
        sample_new_pose.apply(seed);
        
        //check if seed clashes
        AtomPointerVector transformed_seed_atoms = seed->getAtoms();
        seed_clash = isClash(*target_PS, target_BB_atoms, transformed_seed_atoms); //from structgen
        if (attempts > max_attempts) seed_clash = false;
        if (seed_clash) {
            //reset transformations
            sample_new_pose.makeIdentity();
            if (orientation) {
                rotate_around_x.makeIdentity(); rotate_around_y.makeIdentity(); rotate_around_z.makeIdentity();
            }
            if (position) {
                translate_to_bounding_box_origin.makeIdentity(); translate_into_box.makeIdentity();
            } else {
                translate_to_centroid.makeIdentity();
            }
        }
    }
    for (Atom* A_copy : original_seed_atoms_copy) delete A_copy; //the atom copies are no longer necessary
    return attempts;
}


/* --------- naiveSeedsfromDB --------- */

void naiveSeedsFromDB::newPose(string output_path, string out_name, bool position, bool orientation, vector<Residue*> binding_site) {
    
    if ((position == false) && (orientation == false)) MstUtils::error("Neither position or orientation have been selected to be randomized. Nothing to do!");
    cout << "Will randomize: ";
    if (position) cout << "position\t";
    if (orientation) cout << "orientation";
    cout << endl;
    
    //construct bounding box
    structureBoundingBox bounding_box = structureBoundingBox(peptide);
    
    //open new seed binary file
    string seedBinaryPath_out = output_path + out_name + ".bin";
    StructuresBinaryFile seeds_out(seedBinaryPath_out,false);
    
    //open seed data file
    string seedRetry_out = output_path + out_name + "_newpose_attempts.out";
    fstream retry_out;
    MstUtils::openFile(retry_out, seedRetry_out, fstream::out, "naiveSeedsFromBin::newPose");
    //header
    retry_out << "name\tattempts" << endl;
    
    //open seed secondary structure file
    string seedSecStruct_out = output_path + out_name + "_secondary_structure.out";
    fstream secstruct_out;
    MstUtils::openFile(secstruct_out, seedSecStruct_out, fstream::out, "naiveSeedsFromBin::newPose");
    //header
    secstruct_out << "name\tsecondary_structure" << endl;
    
    //randomize seeds and write to file
    int count = 0;
    while (seeds.hasNext()) {
        //read original seed from binary file
        Structure* extended_fragment = seeds.next();
        Chain* seed_C = extended_fragment->getChainByID(seed_chain_id);
        CartesianPoint seed_original_centroid = AtomPointerVector(seed_C->getAtoms()).getGeometricCenter();
        int seed_length = seed_C->residueSize();
        
        //sample new seed from DB
        pair<Structure*,string> seed_data = seedSampler.getSeed(seed_length);
        Structure* new_seed = seed_data.first;
        
        //randomize position/orientation
        int attempts = transform(new_seed, bounding_box, position, orientation, seed_original_centroid);
        
        //write transformed seed to binary file (including meta-data)
        //structure
        seeds_out.appendStructure(new_seed);
        //meta-data
        for (string prop_name : int_properties) {
            int value = seeds.getStructurePropertyInt(prop_name,extended_fragment->getName());
            seeds_out.appendStructurePropertyInt(prop_name,value);
        }
        for (string prop_name : real_properties) {
            mstreal value = seeds.getStructurePropertyReal(prop_name,extended_fragment->getName());
            seeds_out.appendStructurePropertyReal(prop_name,value);
        }
        
        //write attempts to info file
        retry_out << new_seed->getName() << "\t" << attempts << endl;
        
        //write secondary structure
        secstruct_out << new_seed->getName() << "\t" << seed_data.second << endl;
        
        delete extended_fragment;
        delete new_seed;
        count++;
    }
    seeds.reset();
    int num_seeds = seeds.structureCount();
    cout << "There are " << num_seeds << " seeds in the input binary file. After randomization, there are " << count << " seeds in the output binary file" << endl;
    retry_out.close();
    secstruct_out.close();
};

/* --------- seedStatistics --------- */

seedStatistics::seedStatistics(Structure& S, string p_id) : complex(S) {
    neigborhood = 10.0;
    
    //get target structure
    peptide = complex.getChainByID(p_id);
    vector<Residue*> target_residues;
    for (Residue* R : complex.getResidues()) if (R->getChainID() != p_id) target_residues.push_back(R);
    target = Structure(target_residues);
    
    //generate proximity search of target backbone
    cout << "Extracting target backbone to construct a promixity search object..." << endl;
    if (!RotamerLibrary::hasFullBackbone(target)) cout << "warning: target structure is missing backbone atoms!" << endl;
    target_BB_atoms = RotamerLibrary::getBackbone(target);
    target_BB_structure = Structure(target_BB_atoms);
    target_PS = new ProximitySearch(target_BB_atoms, neigborhood/2);
}

void seedStatistics::writeStatisticstoFile(string seedBinaryPath_in, string output_path, string output_name, int num_final_seeds) {
    
    //open seed info file
    string seed_info = output_path + output_name + "_statistics.info";
    fstream seed_out;
    MstUtils::openFile(seed_out, seed_info, fstream::out, "seedStatistics::writeStatisticstoFile");
    //header
    seed_out << "name\tbounding_sphere_radius\tmin_distance_centroid_protein\tseed_length" << endl;
    
    StructuresBinaryFile bin_file(seedBinaryPath_in,true);
    
    long num_seeds = bin_file.structureCount();
    mstreal skip_probability = 1 - min(mstreal(1),mstreal(num_final_seeds)/num_seeds);
    cout << "There are " << num_seeds << " seeds in the input binary file. Skip probability is " << skip_probability << endl;
    
    
    int count = 0;
    while (bin_file.hasNext() == true) {
        mstreal sampled_value = MstUtils::randUnit();
        if (sampled_value <= skip_probability) {
            bin_file.skip();
            continue;
        }
        Structure* extfrag = bin_file.next();
        Chain* seed_C = extfrag->getChainByID("0");
        Structure* seed = new Structure(*seed_C);
        mstreal bounding_sphere_radius = boundingSphereRadius(seed);
        mstreal distance = centroid2NearestProteinAtom(seed);
        seed_out << extfrag->getName() << "\t" << bounding_sphere_radius << "\t" << distance << "\t" << seed->residueSize() << endl;
        delete extfrag;
        delete seed;
        count++;
    }
    cout << "In the end, wrote the info for " << count << " seeds" << endl;
    seed_out.close();
}

mstreal seedStatistics::boundingSphereRadius(Structure *seed) {
    vector<Atom*> seed_atoms = seed->getAtoms();
    CartesianPoint centroid = AtomPointerVector(seed_atoms).getGeometricCenter();
    mstreal max_distance = 0;
    for (Atom* A : seed_atoms) {
        mstreal new_distance = centroid.distance(A);
        if (new_distance > max_distance) max_distance = new_distance;
    }
    return max_distance;
}

mstreal seedStatistics::centroid2NearestProteinAtom(Structure* seed) {
    vector<Atom*> seed_atoms = seed->getAtoms();
    CartesianPoint centroid = AtomPointerVector(seed_atoms).getGeometricCenter();
    return point2NearestProteinAtom(centroid);
}

mstreal seedStatistics::point2NearestProteinAtom(CartesianPoint point) {
    //get nearby points
    vector<int> closest_points;
    mstreal max_d = neigborhood;
    while (closest_points.empty()) {
        closest_points = target_PS->getPointsWithin(point,0.0,max_d);
        max_d = max_d*2;
        if (max_d > 80.0) MstUtils::error("Dude this seed is nowhere near the protein...");
    }
    
    //find the closest one;
    mstreal min_distance = INFINITY;
    for (int i : closest_points) {
        Atom* A = target_BB_atoms[i];
        mstreal new_distance = point.distance(A);
        if (new_distance < min_distance) min_distance = new_distance;
    }
    return min_distance;
}
