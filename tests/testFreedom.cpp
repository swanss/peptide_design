#include "mstsystem.h"
#include "mstoptions.h"
#include "mstcondeg.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("For visualizing the results of the freedom calculation in ConFind.");
    op.addOption("p", "The structure which will be used to compute freedom",true);
    op.setOptions(argc,argv);
    
    string RL = "/Users/sebastianswanson/Keating/utilities/repos/MST/testfiles/rotlib.bin";
    //  string RL = "/home/ifs-users/swans/MST_workspace/MST/testfiles/rotlib.bin";
    
    Structure S(op.getString("p"));
    ConFind C(RL,S,false);
    ConFind C_strict(RL,S,true);
    
    vector<Residue*> all_residues = S.getResidues();
    
    fstream info_out;
    MstUtils::openFile(info_out, "./confind.tsv",fstream::out);
    info_out << "chain_id\tres_num\tfreedom\tcrowdedness\tstrict_freedom\tstrict_crowdedness" << endl;
    
    for (Residue* R : all_residues) {
        mstreal freedom = C.getFreedom(R);
        mstreal crowd = C.getCrowdedness(R);
        mstreal strict_freedom = C_strict.getFreedom(R);
        mstreal strict_crowd = C_strict.getCrowdedness(R);
        cout << *R << " freedom: " << freedom << endl;
        for (Atom* A : R->getAtoms()) A->setB(freedom);
        
        info_out << R->getChainID() << "\t" << R->getNum() << "\t";
        info_out << freedom << "\t";
        info_out << crowd << "\t";
        info_out << strict_freedom << "\t";
        info_out << strict_crowd << endl;
    }
    
    info_out.close();
    
    S.writePDB("test_freedom.pdb");
}

