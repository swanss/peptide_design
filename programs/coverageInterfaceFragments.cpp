//
//  searchInterfaceFragments.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 7/26/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "benchmarkutilities.h"
#include "coverage.h"
#include "secondarystructure.h"
#include "termextension.h"
#include "utilities.h"


int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Finds contacts between the peptide and the protein, uses these to define fragments and looks at the coverage by the peptide aligning segments.");
    op.addOption("pdb", "path to a .pdb file. Can include just protein chains, or both protein/peptide chains", true);
    op.addOption("peptide", "peptide chain ID. Only necessary if the provided .pdb file contains a peptide chain");
    op.addOption("max_rmsd", "The max RMSD threshold used when determining whether a seed aligns to the peptide or not. (default 2.0)");
    op.addOption("config","Path to the configuration file (specifies fasst database and rotamer library)",true);
    op.addOption("flanking_res","The number of residues flanking a contact to include when creating a fragment (default 2).");
    op.setOptions(argc, argv);
    
    MstTimer timer;
    
    // Variables provided by user
    Structure S(op.getString("pdb")), complex;
    string p_cid = op.getString("peptide","");
    mstreal max_rmsd = op.getReal("max_rmsd",2.0);
    configFile config(op.getString("config"));
    int flanking_res = op.getInt("flanking_res",2);
    
    //extract backbone
    RotamerLibrary::extractProtein(complex, S);
    complex.setName(S.getName());
    
    // Make directories
    bool makeParents = true;
    // Fragment output folder
    string outDir = "output/";
    string covDir = "coverage/";
    MstSys::cmkdir(outDir,makeParents);
    MstSys::cmkdir(covDir,makeParents);
    
    selector sel(complex);
    // Set the sequence of the peptide to "unknown"
    vector<Residue*> peptide_res = sel.selectRes("chain "+p_cid);
    cout << "Selected " << peptide_res.size() << " peptide residues to be renamed to 'UNK'" << endl;
    for (Residue* R : peptide_res) {
        R->setName("UNK");
    }
    //Initialize the interface coverage class now, so that the contacts can be used to define interface fragments
    interfaceCoverage IC(complex, p_cid, config.getRL());
    set<pair<Residue*,Residue*>> contactResidues = IC.getContactingResidues();
    
    cout << contactResidues.size() << " unique contacts between the peptide and protein" << endl;
    if (contactResidues.empty()) MstUtils::error("No residues selected for TERM Extension");
    
    cout << "Generate interface fragments and search for matches" << endl;
    
    searchInterfaceFragments SIF(contactResidues,config.getDB());
    SIF.setFlankingResidues(flanking_res);
    
    //Find matches to interface fragments and write to binary file
    SIF.findMatches(outDir);
    
    SIF.writeFragments(outDir);
    
    //Map the seed coverage
    /*
     Vary parameters that filter the seeds
     - match RMSD. The RMSD of the match that generated the seed (depends on the complexity of the fragment)
     - sequence constraint. When applied, only matches with the same residue at the central position are accepted
     */
    
    string interface_matches = "interface_matches";
    string interface_matches_realigned = "interface_matches_repositioned";

    IC.setMaxRMSD(max_rmsd);
    
    IC.writePeptideResidues(covDir + interface_matches + "_");
    IC.writeContacts(covDir + interface_matches + "_");
    
    cout << "Search for segments of seed chains (interface fragment matches) that map to the peptide..." << endl;
    IC.setSeeds(outDir + interface_matches + ".bin");
    IC.findCoveringSeeds();
    cout << "Write coverage to files..." << endl;

    IC.writeAllAlignedSeedsInfo(covDir + interface_matches + "_");
    IC.writeBestAlignedSeeds(covDir + interface_matches + "_",1);
    
    cout << "Search for segments of seed chains (interface fragment matches realigned) that map to the peptide..." << endl;
    IC.setSeeds(outDir + interface_matches_realigned + ".bin");
    IC.findCoveringSeeds();
    cout << "Write coverage to files..." << endl;
    IC.writeAllAlignedSeedsInfo(covDir + interface_matches_realigned + "_");
    IC.writeBestAlignedSeeds(covDir + interface_matches_realigned + "_",1);
    
    cout << "done" << endl;
    return 0;
}
