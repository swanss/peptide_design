#include "freesasaext.h"

void sasaCalculator::setAllowedAA() {
    allowedAA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"};
}

void sasaCalculator::setStructure(Structure* _parentStructure) {
    sasaDefined = false;
    freesasa_structure_free(structure);
    freesasa_node_free(root_node);
    
    parentStructure = _parentStructure;
    structure = freesasa_structure_new();

    // add atoms to freesasa_structure
    vector<Residue*> residues = parentStructure->getResidues();
    for (Residue* R: residues) {
        // only allowedAA have their atom radii defined
        if (allowedAA.find(R->getName()) == allowedAA.end()) {
            cout << "Warning: excluding " << R->getChainID() << R->getNum() << " with residue name " << R->getName() << endl;
            continue;
        }
        for (Atom* A : R->getAtoms()) {
                // exclude hydrogens
                if (RotamerLibrary::isHydrogen(A)) continue;
                
                //not sure how to convert string to char array, so I tried this
                char* atom_name = new char [4];
                string atom_name_str = A->getName();
                leftPadTo(atom_name_str,4);
                strcpy(atom_name,atom_name_str.c_str());

                char* residue_name = new char [3];
                string residue_name_str = A->getResidue()->getName();
                leftPadTo(residue_name_str,3);
                strcpy(residue_name,residue_name_str.c_str());

                char* residue_number = new char [4];
                string residue_number_str = MstUtils::toString(A->getResidue()->getNum());
                leftPadTo(residue_number_str,4);
                strcpy(residue_number,residue_number_str.c_str());

                if (A->getChain()->getID().size() > 1) MstUtils::error("Structure has a chain with an ID longer than a single character");
                string chainID = A->getChain()->getID();
                char chain_label = chainID.at(0);
                
        //        cout << atom_name << " " << residue_name << " " << residue_number << " " << chain_label << endl;
                freesasa_structure_add_atom_wopt(structure,atom_name,residue_name,residue_number,chain_label,A->getX(),A->getY(),A->getZ(),classifier,FREESASA_SKIP_UNKNOWN);
            }
    }
    
    // build a map of residues in the structure
    for (Residue* R : parentStructure->getResidues()) {
        pair<string,int> key(R->getChainID(),R->getNum());
        if (resMap.count(key) > 0) MstUtils::error("Error duplicate position in structure at "+R->getChainID()+","+MstUtils::toString(R->getNum())+" cannot create a map that uniquely identifies each residue in structure","sasaCalculator::setStructure");
        resMap[pair<string,int>(R->getChainID(),R->getNum())] = R;
    }
}

void sasaCalculator::computeSASA() {
    //compute SASA with default parameters and assign output to freesasa_node object
    string name = parentStructure->getName();
    char* name_char = new char[name.size() + 1];
    strcpy(name_char,name.c_str());
    root_node = freesasa_calc_tree(structure,&params,name_char);
}

map<Residue*,mstreal> sasaCalculator::getResidueSASA(bool relative) {
    map<Residue*,mstreal> vals;
    for (Residue* R : parentStructure->getResidues()) vals[R] = 0.0; //default value

    //there should only be one structure in the tree
    freesasa_node *structure_node = freesasa_node_children(freesasa_node_children(root_node));

    // iterate over chains and residues
    freesasa_node* chain_node = freesasa_node_children(structure_node);
    freesasa_node* residue_node = NULL;
    string chainID;
    while (chain_node != NULL) {
        chainID = freesasa_node_name(chain_node);
        residue_node = freesasa_node_children(chain_node);
        string resNum;
        while (residue_node != NULL) {
            resNum = freesasa_node_residue_number(residue_node);
            Residue* R = resMap.at(pair<string,int>(chainID,MstUtils::toInt(resNum)));
            const freesasa_nodearea *residue_area = freesasa_node_area(residue_node);
            if (relative) {
                const freesasa_nodearea *residue_ref_area = freesasa_node_residue_reference(residue_node);
//                cout << chainID << resNum << " " << residue_area->total << " " << residue_ref_area->total << endl;
                vals[R] = mstreal(residue_area->total) / mstreal(residue_ref_area->total);
            } else {
                vals[R] = residue_area->total;
            }
            residue_node = freesasa_node_next(residue_node);
        }
        chain_node = freesasa_node_next(chain_node);
    }

    return vals;
}

void sasaCalculator::leftPadTo(string& str, const size_t num, const char paddingChar) {
    if (str.size() < num) str.insert(0, num - str.size(), paddingChar);
    else if (str.size() > num) MstUtils::error("Provided string value: "+str+" is too large to be padded into a char of size: "+MstUtils::toString(num),"sasaCalculator::padTo");
}
