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
    op.setTitle("Compares the backbone of structure A and B. Computes the RMSD, residue vector angle, and computes contacts.");
    op.addOption("structureA","The first structure that will be considered",true);
    op.addOption("selA","The selection of chains within the first structure",false);
    op.addOption("structureB","The second structure that will be considered",true);
    op.addOption("selB","The selection of chains within the second structure",false);
    op.addOption("resMap","A tsv that specifies the mapping between residues in structure A and B. If provided, residue selections are not used",false);
    op.addOption("config","Path to the configuration file (specifies fasst database and rotamer library)",true);
    op.setOptions(argc, argv);
    
    Structure structureA(op.getString("structureA"));
    Structure structureB(op.getString("structureB"));
    configFile config(op.getString("config"));
    
    string structureAName = MstSys::splitPath(structureA.getName(),1);
    string structureBName = MstSys::splitPath(structureB.getName(),1);
    string fileNamePrefix = structureAName + "_vs_" + structureBName + "_";
    
    Structure selectedStructureA, selectedStructureB;
    set<string> selectedChainsA, selectedChainsB;
    vector<Residue*> selectedResA, selectedResB;
    
    // remove sidechains
    if (!RotamerLibrary::hasFullBackbone(structureA) || !RotamerLibrary::hasFullBackbone(structureB)) {
        MstUtils::error("Structure A or B is missing backbone atoms");
    }
    structureA = Structure(RotamerLibrary::getBackbone(structureA));
    structureB = Structure(RotamerLibrary::getBackbone(structureB));
    
    if (op.isGiven("resMap")) {
        // read res map in, construct vector of chain id/res num for getting residues
        coverageBenchmarkUtils::getResiduesFromMap(op.getString("resMap"),structureA,structureB,selectedResA,selectedResB);
    } else if (op.isGiven("selA") || op.isGiven("selB")) {
        selector selectA(structureA);
        selector selectB(structureB);
                
        selectedResA = selectA.selectRes(op.getString("selA"));
        selectedResB = selectB.selectRes(op.getString("selB"));
    }
    
    selectedStructureA = Structure(selectedResA);
    selectedStructureB = Structure(selectedResB);
    
    cout << "res size A: " << selectedStructureA.residueSize() << " res size B: " << selectedStructureB.residueSize() << endl;
    cout << "chain size A: " << selectedStructureA.chainSize() << " chain size B: " << selectedStructureB.chainSize() << endl;
    for (int i = 0; i < selectedStructureA.chainSize(); i++) selectedChainsA.insert(selectedStructureA.getChain(i).getID());
    for (int i = 0; i < selectedStructureB.chainSize(); i++) selectedChainsB.insert(selectedStructureB.getChain(i).getID());
        
    // check that everything worked
    MstUtils::assertCond((selectedStructureA.atomSize() == selectedStructureB.atomSize()),"Different number of atoms in selected structure: "+MstUtils::toString(selectedStructureA.atomSize())+" and "+MstUtils::toString(selectedStructureB.atomSize()));
    MstUtils::assertCond((selectedStructureA.residueSize() == selectedStructureB.residueSize()),"Different number of residues in selected structure: "+MstUtils::toString(selectedStructureA.residueSize())+" and "+MstUtils::toString(selectedStructureB.residueSize()));
    MstUtils::assertCond(((selectedChainsA.size() > 0) || (selectedChainsB.size() > 0)),"Missing chains in selected structure: "+MstUtils::toString(selectedChainsA.size())+" and "+MstUtils::toString(selectedChainsB.size()));
    
    // now compute the statistics
    mstreal rmsd = coverageBenchmarkUtils::writeRMSDtoFile(fileNamePrefix,selectedStructureB,selectedStructureA);
    cout << "The total RMSD between the structure A and B is: " << rmsd << endl;
    
    // get the contacts for both structures (doesn't matter which chains we select as peptide
    interfaceCoverage ICA(structureA, *selectedChainsA.begin(), config.getRL());
    interfaceCoverage ICB(structureB, *selectedChainsB.begin(), config.getRL());

    coverageBenchmarkUtils::writeContactstoFile(structureAName+"_", &ICA, structureA, selectedChainsA, config.getRL());
    coverageBenchmarkUtils::writeContactstoFile(structureBName+"_", &ICB, structureB, selectedChainsB, config.getRL());
    
    return 0;
}
