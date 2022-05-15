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
    op.addOption("getTarget", "The target index of a structure you want to save as PDB");
    op.addOption("mem_save","0: no memory save, 1: no sidechains, 2: destroy structures after getting atoms (default: 0)",false);
    op.addOption("stride","If provided, will write all STRIDE secondary classifications to a csv file",false);
    op.setOptions(argc, argv);
    
    int mem_save = op.getInt("mem_save",0);
    bool stride = op.isGiven("stride");
    int targetID2get = op.getInt("getTarget",-1);
    
    //import DB
    FASST F;
    cout << "Reading DB..." << endl;
    F.readDatabase(op.getString("dbPath"),mem_save);
    cout << "Done reading DB" << endl;
    
    //Structure properties
    int num_target = F.numTargets();
    int num_chain = 0;
    int num_res = 0;
    for (int target_id = 0; target_id < num_target; target_id++) {
        Structure* S = F.getTarget(target_id);
        num_chain += S->chainSize();
        num_res += S->residueSize();
    }
    
    cout << "dbPath: " << op.getString("dbPath") << endl;
    cout << "Structures in DB: " << num_target << endl;
    cout << "Chains in DB: " << num_chain << endl;
    cout << "Residues in DB: " << num_res << endl;
    
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
        
        out << "target_idx,target_name,res_idx,chain_id,res_num,stride,res_name" << endl;
        
        for (int i = 0; i < F.numTargets(); i++) {
            Structure* target = F.getTarget(i);
            string target_name = MstSys::splitPath(target->getName(),1);
            for (Residue* R : target->getResidues()) {
                int idx = R->getResidueIndex();
                out << i << ",";
                out << target_name << "," << idx << "," << R->getChainID() << ",";
                out << R->getNum() << "," << F.getResidueStringProperty(i, "stride", idx) << ",";
                out << R->getName() << endl;
            }
        }
        out.close();
    }
    
    if (targetID2get != -1) {
        Structure* target = F.getTarget(targetID2get);
        target->writePDB(MstUtils::toString(targetID2get)+".pdb");
    }
    
    return 0;
}
