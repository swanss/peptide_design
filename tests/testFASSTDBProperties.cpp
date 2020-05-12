//
//  testFASSTDBProperties.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 10/10/19.
//

#include "msttypes.h"
#include "mstfasst.h"
#include "mstoptions.h"

using namespace std;
using namespace MST;

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Reports which properties are defined in a FASSTDB");
    op.addOption("dbPath","The path to the DB",true);
    op.setOptions(argc, argv);
    
    //import DB
    FASST F;
    F.readDatabase(op.getString("dbPath"),1);
    
    cout << "Residue property: phi = " << F.isResiduePropertyDefined("phi") << endl;
    cout << "Residue property: psi = " << F.isResiduePropertyDefined("psi") << endl;
    cout << "Residue property: omega = " << F.isResiduePropertyDefined("omega") << endl;
    cout << "Residue pair property: cont = " << F.isResiduePairPropertyPopulated("cont") << endl;
    cout << "Residue pair property: interfering = " << F.isResiduePairPropertyPopulated("interfering") << endl;
    cout << "Residue pair property: interfered = " << F.isResiduePairPropertyPopulated("interfered") << endl;
    cout << "Residue pair property: bb = " << F.isResiduePairPropertyPopulated("bb") << endl;
    cout << "Residue relationship property sim: " << MstUtils::toString(F.isResidueRelationshipPopulated("sim")) << endl;
    
    return 1;
}
