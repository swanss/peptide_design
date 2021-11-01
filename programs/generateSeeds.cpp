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
    op.addOption("no_seeds","If provided will skip generating seeds");
    op.addOption("only_store_covering","If provided will only write the seeds that are covering, e.g. have some segment aligning to the peptide");
//    op.addOption("adaptive_fragments","Fragments grow as large as possible while still having the required number of matches");
//    op.addOption("disjoint_segments","Fragments grow by adding disjoint segments (if not provided, fragments only grow by adding flanking residues)");
    op.addOption("write_all_files","Writes additional files (helpful for making figures and diagnosing issues)");
    op.setOptions(argc, argv);
        
    MstTimer timer;
    
    // Variables provided by user
    Structure target(op.getString("pdb"));
    string params_file_path = op.getString("params_file");
    string p_cid = op.getString("peptide","");
    string sel_str = op.getString("sel","");
    bool no_seeds = op.isGiven("no_seeds");
    bool only_store_covering = op.isGiven("only_store_covering");
//    bool adaptive_fragments = op.isGiven("adaptive_fragments");
//    bool disjoint_segments = op.isGiven("disjoint_segments");
    bool write_all_files = op.isGiven("write_all_files");
  
    // Open params file
    TEParams params(params_file_path);
    params.printValues();
    configFile config(params.getConfigFile());
    
    // Make directories
    bool makeParents = true;
    // Fragment output folder
    string outDir = "output/";
    MstSys::cmkdir(outDir,makeParents);
    string filesDir = "files/";
    string extFragDir = filesDir + "extFrag/";
    string wholeMatchProteinDir = filesDir + "wholeMatchProtein/";
    string seedDir = filesDir + "seedDir/";
    if (write_all_files) {
        MstSys::cmkdir(extFragDir,makeParents);
        MstSys::cmkdir(wholeMatchProteinDir,makeParents);
        MstSys::cmkdir(seedDir,makeParents);
    }
    
    vector<Residue*> bindingSiteRes;
    interfaceCoverage* IC = nullptr;
    if ((p_cid != "") && (sel_str != "")) {
        // Delete the peptide and then select residues
        Chain* C = target.getChainByID(p_cid);
        target.deleteChain(C);
        selector sel(target);
        bindingSiteRes = sel.selectRes(sel_str);
    } else if (p_cid != "") {
        // Set the sequence of the peptide to "unknown"
        selector sel(target);
        vector<Residue*> peptide_res = sel.selectRes("chain "+p_cid);
        cout << "Selected " << peptide_res.size() << " peptide residues to be renamed to 'UNK'" << endl;
        for (Residue* R : peptide_res) {
            R->setName("UNK");
        }
        //Initialize the interface coverage class now, so that the binding site residues can be used in TERM Extension
        IC = new interfaceCoverage(target, p_cid, config.getRL());
        bindingSiteRes = IC->getBindingSiteRes();
    } else if (sel_str != "") {
        selector sel(target);
        bindingSiteRes = sel.selectRes(sel_str);
    } else {
        // select all residues with a freedom above the threshold
        cout << "Select residues with freedom above the threshold: " << params.freedom_cutoff << endl;
        ConFind C(config.getRL(),target);
        for (Residue* R : target.getResidues()) {
            mstreal freedom = C.getFreedom(R);
            if (freedom > params.freedom_cutoff) bindingSiteRes.push_back(R);
        }
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
    if (write_all_files) TE.setWriteAllFiles(extFragDir,wholeMatchProteinDir,seedDir);
    timer.start();
    TE.generateFragments();
    timer.stop();
    cout << timer.getDuration() << " seconds to generate fragments" << endl;
    
    timer.start();
    if (only_store_covering) TE.setIC(IC);
    if (!no_seeds) TE.extendFragmentsandWriteStructures(seedTERM::MANY_CONTACT,outDir);
    timer.stop();
    cout << timer.getDuration() << " seconds to generate seeds" << endl;
    
    //write fragments
    cout << "Writing fragment pdbs..." << endl;
    TE.writeFragmentPDBs(outDir);
    TE.writeFragmentClassification(outDir);
    
    if (IC != nullptr) delete IC;
    
    cout << "Done" << endl;
    return 0;
}
