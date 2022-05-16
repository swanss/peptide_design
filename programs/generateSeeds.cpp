//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//peptide_design dependencies
#include "utilities.h"
#include "coverage.h"
#include "termextension.h"
#include "secondarystructure.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Generates segments of protein backbone, or 'interface seeds', around a target protein.");
    op.addOption("targetPDB", "Path to a PDB file of the target protein", true);
    op.addOption("peptideChainID", "Single letter peptide chain ID. Will define a peptide binding site and remove the peptide before generating seeds. Only necessary if the targetPDB contains a peptide chain that should be removed");
    op.addOption("targetSel","A selection string specifying protein residues to generate seeds around. Only will be used if peptideChainID is not provided. If neither are provided, will try to generate seeds around all residues. Ex: 'chain A and resid 122-130'");
    op.addOption("paramsFile","Path to the configuration file (specifies FASST structure database and rotamer library)",true);
    op.addOption("noSeeds","If provided will skip generating seeds");
    op.addOption("onlyStoreCovering","If provided will only write the seeds that are covering, e.g. have some segment aligning to the peptide");
    op.addOption("writeAllFiles","Writes additional files (helpful for making figures and diagnosing issues)");
    op.setOptions(argc, argv);
        
    MstTimer timer;
    
    // Variables provided by user
    Structure target(op.getString("targetPDB"),"ALLOW ILE CD1");
    string params_file_path = op.getString("paramsFile");
    string p_cid = op.getString("peptideChainID","");
    string sel_str = op.getString("targetSel","");
    bool no_seeds = op.isGiven("noSeeds");
    bool only_store_covering = op.isGiven("onlyStoreCovering");
    bool write_all_files = op.isGiven("writeAllFiles");
  
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
    if (only_store_covering) {
        cout << "Will only store seeds that cover the peptide" << endl;
        TE.setIC(IC);
    }
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
