//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "utilities.h"
#include "coverage.h"
#include "benchmarkutilities.h"
//#include "termextension.h"
#include "secondarystructure.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Writes out some basic residue-level info for a protein");
    op.addOption("PDB","The path to the PDB structure",false);
    op.addOption("PDB_list","The path to a file where each line is the path to a PDB structure",false);
    op.setOptions(argc, argv);
    
    string pdb_path = op.getString("PDB","");
    string pdb_list_path = op.getString("PDB_list","");
    
    vector<string> PDB_paths;
    if (!pdb_list_path.empty()) {
        PDB_paths = MstUtils::fileToArray(pdb_list_path);
    } else if (!pdb_path.empty()) {
        PDB_paths.push_back(pdb_path);
    } else {
        MstUtils::error("Must provide --PDB or PDB_list","WriteStructureResidueProperties");
    }
    
    for (string pdb_path : PDB_paths) {
        string pdb = MstSystemExtension::fileName(pdb_path);
        string pdb_name = pdb.substr(0,pdb.length()-4);
        
        cout << "Getting info for structure: " << pdb_name << endl;
        
        Structure S(pdb_path);
        
        fstream out;
        MstUtils::openFile(out, pdb_name+".csv",fstream::out);
        out << "PDB_name,chain_id,res_num,res_idx,res_name,phi,psi" << endl;
        
        for (Residue* R : S.getResidues()) {
            out << pdb_name << "," << R->getChainID() << "," << R->getNum() << ",";
            out << R->getResidueIndex() << "," << R->getName() << ",";
            if (RotamerLibrary::hasFullBackbone(R)) {
                out << R->getPhi() << "," << R->getPsi();
            } else {
                out << 999 << "," << 999;
            }
            out << R->getPhi() << "," << R->getPsi();
            out << endl;
        }
    }

    return 0;
};
