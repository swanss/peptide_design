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
    op.addOption("config","Path to the configuration file (specifies fasst database and rotamer library)",true);
    op.addOption("flanking_res","The number of residues flanking a contact to include when creating a fragment (default 2).");
    op.addOption("match_req", "The fragmenter will attempt to create the largest (ranked by number of residues) fragments that have at least this many matches. During TERM Extension, even if the fragment has more than this number match_num_req matches, only this number will be used to generate seeds. If not defined, defaults to CEN_RES.");
    op.addOption("variable_rmsd", "If match_req == true and this option is provided fragments will not grow, instead the RMSD cutoff will be increased until sufficient matches have be found");
    op.addOption("match_rmsd", "Sets the max rmsd allowed when searching for matches to protein fragments.");
    op.addOption("no_adaptive_rmsd","If provided, will not use adaptive RMSD cutoff.");
    op.addOption("seq", "Constrain the matches to those that have the same central amino acid");
    op.setOptions(argc, argv);
    
    if (op.isGiven("peptide") == op.isGiven("sel")) MstUtils::error("Either a peptide chain ID or a selection string must be provided, but not both");
    
    MstTimer timer;
    
    // Variables provided by user
    Structure target(op.getString("pdb"));
    configFile config(op.getString("config"));
    string p_cid = op.getString("peptide","");
    string sel_str = op.getString("sel","");
    int flanking_res = op.getInt("flanking_res",2);
    mstreal match_rmsd = op.getReal("match_rmsd",1.2);
    bool variable_rmsd = op.isGiven("variable_rmsd");
    bool adaptive_rmsd = !op.isGiven("no_adaptive_rmsd");
    
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
    
    TermExtension TE(config.getDB(), config.getRL(), bindingSiteRes);
    TE.setFlankingNumber(flanking_res);
    TE.setMaxRMSD(match_rmsd);
    TE.setAdaptiveRMSD(adaptive_rmsd);
    if (op.isGiven("seq")) TE.setSeqConst(true);
    timer.start();
    if (op.isGiven("match_req")) {
        cout << "Match requirement: " << op.getInt("match_req") << endl;
        TE.setMatchReq(op.getInt("match_req"));
        if (variable_rmsd) {
            TE.setAdaptiveRMSD(false);
            TE.generateFragments(TermExtension::MATCH_NUM_REQ_CUTOFF);
        } else {
            TE.generateFragments(TermExtension::MATCH_NUM_REQ_SIZE);
        }
    }
    else {
        cout << "No match requirement, CEN_RES mode" << endl;
        TE.generateFragments(TermExtension::CEN_RES);
    }
    timer.stop();
    cout << timer.getDuration() << " seconds to generate fragments" << endl;
    
    timer.start();
    TE.extendFragmentsandWriteStructures(seedTERM::MANY_CONTACT,outDir);
    timer.start();
    cout << timer.getDuration() << " seconds to generate seeds" << endl;
    
    //write fragments
    cout << "Writing fragment pdbs..." << endl;
    TE.writeFragmentPDBs(outDir);
    TE.writeFragmentClassification(outDir);
    
    //write parameters
    TE.storeParameters("parameters.info");
    
    if (IC != nullptr) delete IC;
    
    cout << "done" << endl;
    return 0;
}
