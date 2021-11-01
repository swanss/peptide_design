#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"
#include "msttypes.h"

#include "coverage.h"
#include "freesasaext.h"
#include "utilities.h"
#include "termextension.h"
//#include "scoretargetresidues.h"

void replaceMSE(Structure& S) {
    for (Residue* R : S.getResidues()) {
        if (R->getName() == "MSE") R->setName("MET");
        for (Atom* A : R->getAtoms()) {
            if (A->getName() == "SE") A->setName("SD");
        }
    }
}

int main(int argc, char* argv[]) {
    MstOptions op;
    op.setTitle("Scores the 'bind-ability' of each residue on the target protein.");
    op.addOption("pdb", "Path to a .pdb file. Can include just protein chains, or both protein/peptide chains");
    op.addOption("gz","Path to pdb.gz file. Alternative to 'pdb'");
    op.addOption("paramsFile","Path to the configuration file (specifies fasst database and rotamer library)",true);
    op.addOption("cR", "Nearby seed count radius. Will use this distance to count nearby seeds");
    op.addOption("vR", "Volume calculation radius. Will use this distance to determine the fraction of unoccupied volume");
    op.addOption("expectedNearbyCountsHist","Path to a binned data file. Reports the number of nearby seed residues per match when residues are binned by fraction of unoccupied volume",true);
    op.addOption("expectedGeneratedCountsHist","Path to a binned data file. Reports the number of generated seed residues per match when residues are binned by fraction of unoccupied volume",true);
    op.addOption("protSel","A selection string used to select the residues of the protein");
    op.addOption("peptideChain","The chain ID of the peptide. Required if a peptide-protein complex is provided");
    op.addOption("seeds","Path to a structures binary file containing previously generated seeds. If provided, will skip the seed generating step and count these instead");
    op.setOptions(argc, argv);
    
    MstTimer timer;
    
    string pdb_name = op.getString("pdb","");
    string gz_name = op.getString("gz","");
    string params_file_path = op.getString("paramsFile");
    mstreal cR = op.getReal("cR",12.0);
    mstreal vR = op.getReal("vR",12.0);
    string expectedNearbyCountsHist = op.getString("expectedNearbyCountsHist","");
    string expectedGeneratedCountsHist = op.getString("expectedGeneratedCountsHist","");
    string protSel = op.getString("protSel");
    string peptideChainID = op.getString("peptideChain","");
    string seedBinPath = op.getString("seeds","");
    
    if ((pdb_name == "") && (gz_name == "")) MstUtils::error("Either --pdb or --gz must be provided");
  
    // Open params file
    TEParams params(params_file_path);
    params.printValues();
    configFile config(params.getConfigFile());
    
    // Make directories
    bool makeParents = true;
    string termExt_outDir = "";
    string seedSample_outDir = "";
    if (seedBinPath == "") {
        termExt_outDir = "termext_output/";
        MstSys::cmkdir(termExt_outDir,makeParents);
        seedSample_outDir = "seedsample_output/";
        MstSys::cmkdir(seedSample_outDir,makeParents);
    }
    string resScore_outDir = "scoreres_output/";
    MstSys::cmkdir(resScore_outDir,makeParents);
    
    // Open target pdb, select the binding site residues if applicable, and delete the peptide chain
    /*
     Missing backbone atoms will result in an error, but are reasonably common in structures from
     'the wild' (read: the PDB). To avoid issues we will create a copy of the structure: including
     only those residues with complete backbones. In case there are any issues the PDB will be written
     out so that the user can diagnose if something undesirable has occurred. Note that this structure
     is a "shallow copy", in that all of its atoms point to those in the original structure.
     */
    Structure structure;
    if (pdb_name != "") {
        Structure structure_orig(pdb_name,"ALLOW ILE CD1 QUIET");
        cout << structure_orig.getResidues().size() << " residues in original structure" << endl;
        RotamerLibrary::extractProtein(structure, structure_orig);
    } else {
        // if it is a .gz file, it was likely downloaded from the PDB and needs to be cleaned up
        string tmp = "tmp.pdb";
        MstSys::csystem("gunzip < " + gz_name + " > " + tmp);
        Structure structure_orig(tmp, "ALLOW ILE CD1 QUIET"), structure_extract;
        RotamerLibrary::extractProtein(structure_extract,structure_orig);
        structure_extract.reassignChainsByConnectivity(structure);
        structure.renumber();
    }
    //lazy quick fix
    replaceMSE(structure);
    
    interfaceCoverage* IC = nullptr;
    vector<Residue*> selectedProteinRes;
    if (peptideChainID != "") {
        cout << "Peptide chain provided: '" << peptideChainID << "'. Will use peptide contacts to define a binding site" << endl;
        selector sel(structure);
        vector<Residue*> peptide_res = sel.selectRes("chain "+peptideChainID);
        cout << "Selected " << peptide_res.size() << " peptide residues" << endl;
        //Initialize the interface coverage class now, so that the binding site residues can be used in TERM Extension
        IC = new interfaceCoverage(structure, peptideChainID, config.getRL());
        Structure* target = IC->getTargetStructure();
        if (protSel != "") {
            cout << "Protein selection string provided: " << protSel << " will define the target chains" << endl;
            selector sel(*target);
            selectedProteinRes = sel.selectRes(protSel);
        } else {
            cout << "No protein selection string provided, all non-peptide residues will be selected" << endl;
            selectedProteinRes = target->getResidues();
        }
    } else {
        cout << "No peptide selection chain provided" << endl;
        if (protSel != "") {
            cout << "Protein selection string provided: " << protSel << " will define the target chains" << endl;
            selector sel(structure);
            selectedProteinRes = sel.selectRes(protSel);
        } else {
            cout << "No protein selection string provided, all residues will be selected" << endl;
            selectedProteinRes = structure.getResidues();
        }
    }

    cout << "Selected " << selectedProteinRes.size() << " residues" << endl;
    
    if (seedBinPath == "") {
        cout << "No seeds provided, will generate new seeds around the protein" << endl;
        // Build fragments and use these to generate seeds
        TermExtension TE(config.getDB(), config.getRL(), selectedProteinRes, params);
        TE.setVCalcRadius(vR);
        timer.start();
        TE.generateFragments();
        timer.stop();
        cout << timer.getDuration() << " seconds to generate fragments" << endl;
        
        timer.start();
        TE.extendFragmentsandWriteStructures(seedTERM::MANY_CONTACT,termExt_outDir);
        timer.stop();
        cout << timer.getDuration() << " seconds to generate seeds" << endl;
        
        // Write fragments
        cout << "Writing fragment pdbs..." << endl;
        TE.writeFragmentPDBs(termExt_outDir,expectedGeneratedCountsHist); //not the right path, but will use for now
        TE.writeFragmentClassification(termExt_outDir);
        
        seedBinPath = termExt_outDir + "extendedfragments.bin";
        
        // Sample seeds and write out as PDBs
        mstreal rate = 0.1;
        StructuresBinaryFile seeds(seedBinPath);
        while (seeds.hasNext()) {
            if (rand() / (float)RAND_MAX >= rate) {
                seeds.skip();
                continue;
            }
            Structure *s = seeds.next();
            Structure S_seed;
            extractChains(*s, "0", S_seed);
            S_seed.writePDB(seedSample_outDir+s->getName()+".pdb");
            delete s;
        }
    }

    // Score residues
    cout << "Scoring residues..." << endl;
    set<Residue*> bindingSiteRes;
    
    timer.start();
    nearbySeedScore* scorer;
    if (peptideChainID != "") {
        bindingSiteRes = IC->getBindingSiteResByDistance();
        scorer = new nearbySeedScore(selectedProteinRes, seedBinPath, expectedNearbyCountsHist, params.match_req, cR, vR);
    } else {
        scorer = new nearbySeedScore(selectedProteinRes, seedBinPath, expectedNearbyCountsHist, params.match_req, cR, vR);
    }
    scorer->scoreResidues();
    scorer->writeResidueInfo(resScore_outDir,bindingSiteRes);
    timer.stop();
    cout << timer.getDuration() << " seconds to score residues and write out results" << endl;
    delete scorer;
    delete IC;
    
    cout << "Done" << endl;
    return 0;
}
