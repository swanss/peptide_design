#include "msttypes.h"
#include "mstfasst.h"
#include "mstoptions.h"
#include "mstsystem.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Reports which properties are defined in a FASSTDB");
    op.addOption("dbPath","The path to the DB",true);
    op.addOption("mem_save","0: no memory save, 1: no sidechains, 2: destroy structures after getting atoms (default: 0)",false);
    op.addOption("stride","If provided, will write all STRIDE secondary classifications to a csv file",false);
    op.setOptions(argc, argv);
    
    int mem_save = op.getInt("mem_save",0);
    bool stride = op.isGiven("stride");
    
    //import DB
    FASST F;
    F.readDatabase(op.getString("dbPath"),mem_save);
    
    //Structure properties
    int target_num = F.numTargets();
    int chain_num = 0;
    for (int target_id = 0; target_id < target_num; target_id++) {
        Structure* S = F.getTarget(target_id);
        chain_num += S->chainSize();
    }
    
    cout << "dbPath: " << op.getString("dbPath") << endl;
    cout << "Structures in DB: " << target_num << endl;
    cout << "Chains in DB: " << chain_num << endl;
    
    cout << "Residue property: phi = " << F.isResiduePropertyDefined("phi") << endl;
    cout << "Residue property: psi = " << F.isResiduePropertyDefined("psi") << endl;
    cout << "Residue property: omega = " << F.isResiduePropertyDefined("omega") << endl;
    cout << "Residue property: env = " << F.isResiduePropertyDefined("env") << endl;
    cout << "Residue property: stride = " << F.isResidueStringPropertyDefined("stride") << endl;
    cout << "Residue pair property: cont = " << F.isResiduePairPropertyPopulated("cont") << endl;
    cout << "Residue pair property: interfering = " << F.isResiduePairPropertyPopulated("interfering") << endl;
    cout << "Residue pair property: interfered = " << F.isResiduePairPropertyPopulated("interfered") << endl;
    cout << "Residue pair property: bb = " << F.isResiduePairPropertyPopulated("bb") << endl;
    cout << "Residue relationship property sim: " << MstUtils::toString(F.isResidueRelationshipPopulated("sim")) << endl;
    
    // create map of prop names for contact degree with amino acid constraints
    vector<string> aaNames = {"ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
        "LYS", "MET", "PHE", "SER", "THR", "TRP", "TYR", "VAL", "ALA"}; //aa allowed in rotamer library
    map<string, string> aaToProp;
    for (string aa : aaNames) aaToProp[aa] = "cont"+aa;
    
    for (auto it: aaToProp) cout << "Residue pair property: " << it.second << " = " << F.isResiduePairPropertyPopulated(it.second) << endl;
    
    if (stride) {
        cout << "Writing out stride info..." << endl;
        
        fstream out;
        string output_path = "DB_STRIDE_info.csv";
        MstUtils::openFile(out, output_path, fstream::out);
        
        out << "target,res_idx,chain_id,res_num,stride" << endl;
        
        for (int i = 0; i < F.numTargets(); i++) {
            Structure* target = F.getTarget(i);
            string target_name = MstSys::splitPath(target->getName(),1);
            for (Residue* R : target->getResidues()) {
                int idx = R->getResidueIndex();
                out << target_name << "," << idx << "," << R->getChainID() << ",";
                out << R->getNum() << "," << R->getNum() << "," << F.getResidueStringProperty(i, "stride", idx);
                out << endl;
            }
        }
        
        out.close();
    }
    
    return 0;
}
