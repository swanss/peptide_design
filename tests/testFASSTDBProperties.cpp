 #include "msttypes.h"
#include "mstfasst.h"
#include "mstoptions.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Reports which properties are defined in a FASSTDB");
    op.addOption("dbPath","The path to the DB",true);
    op.addOption("secStruct","If given, assigns a secondary classification to each residue of each structure in the database",true);
    op.addOption("mem_save","0: no memory save, 1: no sidechains, 2: destroy structures after getting atoms",true);
    op.setOptions(argc, argv);
    
    int mem_save = op.getInt("mem_save");
    
    //import DB
    FASST F;
    F.readDatabase(op.getString("dbPath"),mem_save);
    
    if (op.isGiven("secStruct")) {
        
    }
    
    //Structure properties
    int target_num = F.numTargets();
    int chain_num = 0;
    for (int target_id = 0; target_id < target_num; target_id++) {
        Structure* S = F.getTarget(target_id);
        chain_num += S->chainSize();
    }

    vector<string> allProperties = F.getAllResidueProperties();
    cout << "All Properties" << endl;

    for (string property : allProperties) {
        cout << property << endl;
    }

    cout << endl << endl;
    
    cout << "dbPath: " << op.getString("dbPath") << endl;
    cout << "Structures in DB: " << target_num << endl;
    cout << "Chains in DB: " << chain_num << endl;
    
    cout << "Residue property: phi = " << F.isResiduePropertyDefined("phi") << endl;
    cout << "Residue property: psi = " << F.isResiduePropertyDefined("psi") << endl;
    cout << "Residue property: omega = " << F.isResiduePropertyDefined("omega") << endl;
    cout << "Residue property: env = " << F.isResiduePropertyDefined("env") << endl;
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
    
    return 0;
}
