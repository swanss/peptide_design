//
// TERMExtensionBenchmark.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 2/24/19.
//


//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "utilities.h"
#include "coverage.h"
#include "termextension.h"
#include "secondarystructure.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Generates segments of protein backbone, referred to as 'seeds', that are potentially compatible with the surface of a target protein through TERM Extension. These are saved in a binary file.");
    op.addOption("pdb", "path to a .pdb file. Can include just protein chains, or both protein/peptide chains", true);
    op.addOption("peptide", "peptide chain ID. Only necessary if the provided .pdb file contains a peptide chain");
    op.addOption("sel","a selection string that specifies the protein residues to generate seeds around. Necessary if the provided PDB file does not include peptide chains");
    op.addOption("params_file","Path to the configuration file (specifies fasst database and rotamer library)",true);
    op.addOption("only_store_covering","If provided, will only write the seeds that are covering, e.g. have some segment aligning to the peptide");
    op.addOption("adaptive_fragments","Fragments grow as large as possible while still having the required number of matches");
    op.setOptions(argc, argv);
    
    if (op.isGiven("peptide") == op.isGiven("sel")) MstUtils::error("Either a peptide chain ID or a selection string must be provided, but not both");
    
    MstTimer timer;
    
    // Variables provided by user
    Structure target(op.getString("pdb"));
    string params_file_path = op.getString("params_file");
    string p_cid = op.getString("peptide","");
    string sel_str = op.getString("sel","");
    bool only_store_covering = op.isGiven("only_store_covering");
    bool adaptive_fragments = op.isGiven("adaptive_fragments");
  
    // Open params file
    TEParams params(params_file_path);
    params.printValues();
    configFile config(params.getConfigFile());
    
    // Make directories
    bool makeParents = true;
    // Fragment output folder
    string outDir = "output/";
    MstSys::cmkdir(outDir,makeParents);
    
    selector sel(target);
    vector<Residue*> bindingSiteRes;
    interfaceCoverage* IC = nullptr;
    if (p_cid != "") {
        // Set the sequence of the peptide to "unknown"
        vector<Residue*> peptide_res = sel.selectRes("chain "+p_cid);
        cout << "Selected " << peptide_res.size() << " peptide residues to be renamed to 'UNK'" << endl;
        for (Residue* R : peptide_res) {
            R->setName("UNK");
        }
        //Initialize the interface coverage class now, so that the binding site residues can be used in TERM Extension
        IC = new interfaceCoverage(target, p_cid, config.getRL());
        bindingSiteRes = IC->getBindingSiteRes();
    } else {
        bindingSiteRes = sel.selectRes(sel_str);
    }
    
    if (bindingSiteRes.empty()) MstUtils::error("No residues selected for TERM Extension");
    cout << "The following residues were selected for TERM Extension: ";
    for (Residue* R : bindingSiteRes) cout << R->getChainID() << R->getNum() << " ";
    cout << endl;
    
    // Write out the secondary structure of each residue in the original protein
    fstream out;
    string output_path = "secondarystructure.tsv";
    MstUtils::openFile(out, output_path, fstream::out);
    
    secondaryStructureClassifier classifier;
    classifier.writeResClassificationtoTSV(target, out);
    out.close();
    
    TermExtension TE(config.getDB(), config.getRL(), bindingSiteRes, params);
    timer.start();
    if (adaptive_fragments) TE.generateFragments(TermExtension::ADAPTIVE_SIZE);
    else TE.generateFragments(TermExtension::CEN_RES);
    timer.stop();
    cout << timer.getDuration() << " seconds to generate fragments" << endl;
    
    timer.start();
    if (only_store_covering) TE.setIC(IC);
    TE.extendFragmentsandWriteStructures(seedTERM::MANY_CONTACT,outDir);
    timer.stop();
    cout << timer.getDuration() << " seconds to generate seeds" << endl;
    
    //write fragments
    cout << "Writing fragment pdbs..." << endl;
    TE.writeFragmentPDBs(outDir);
    TE.writeFragmentClassification(outDir);
    
    //write parameters
    TE.storeParameters("parameters.info");
    
    if (IC != nullptr) delete IC;
    
    cout << "Done" << endl;
    return 0;
}
