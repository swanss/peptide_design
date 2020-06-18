//
//  testTERMExtensionVarySequence.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 3/2/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "utilities.h" //for generateAllCombinationsKRes
#include "coverage.h" //to directly assess coverage
#include "termextension.h"
#include "secondarystructure.h"
#include "structure_iter.h"


int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Takes a single fragment from a structure and expands a selected residue by TERM Extension. This is repeated, 20 times, with a different amino acid placed at the selected residue as a search constraint");
  op.addOption("pdb", "Structure file in which the TERM exists", true);
  op.addOption("sel", "The central esidue to expand in the TERM", true);
  op.addOption("match_num_req", "The fragmenter will attempt to create the largest (containing the most residues) fragments that have at least this many matches. During TERM Extension, even if the fragment has more than this number match_num_req matches, only this number will be used to generate seeds.  If not defined, defaults to CEN_RES.");
  op.addOption("complexity_scan", "If provided uses COMPLEXITY_SCAN option for making fragments, otherwise defaults to CEN_RES");
  op.addOption("config", "Path to configuration file (specifies fasst database and rotamer library)",true);
  op.setOptions(argc, argv);
  
  configFile config(op.getString("config"));
  
  int match_num_req = op.getInt("match_num_req");

  
  // Import structure
  Structure target(op.getString("pdb"));
  
  // Select peptide residues, use to get all peptide chain pointers
  selector sel(target);
  vector<Residue*> res = sel.selectRes(op.getString("sel"));
  if (res.size() != 1) MstUtils::error("Only one residue should be selected from the target protein, " + MstUtils::toString(res.size()) + " currently selected.");
  cout << "Selected " << *res.front() << endl;
  
  //// Make directories
  // Extended fragment output folder
  bool makeParents = true;
  string extfragDir = "termextension_output/";
  MstSys::cmkdir(extfragDir,makeParents);
  
  MstTimer timer;
  vector<Structure*> extended_fragments;
  
  //// initialize fragmenter
  cout << "Initializing fragmenter object..." << endl;
  TermExtension F(config.getDB(), config.getRL(), res);
  
  //Store the settings of the fragmenter/contacts
  F.setSeqConst(true);
  F.storeParameters("fragmenter.info");
  
  Residue* cen_res = res.front();
  //// Try resetting the central position to all natural amino acids and generating extended fragments
  SeqTools seqtool;
  vector<string> aa3 = seqtool.getAA3();
  
  //// Open file
  string aa_info = "aa.tsv";
  fstream fs;
  MstUtils::openFile(fs,aa_info,fstream::out);
  fs << "amino_acid\tnum_matches" << endl;
  
  secondaryStructureClassifier classifier;
  for (string aa : aa3) {
    cen_res->setName(aa);
    
    ////Create fragments
    timer.start();
    if (op.isGiven("match_num_req")) {
      F.setMatchReq(match_num_req);
      cout << "Generating fragments in the match number requirement mode..." << endl;
      F.generateFragments(TermExtension::MATCH_NUM_REQ);
    } else if (op.isGiven("complexity_scan")) {
      cout << "Generating fragments in the complexity scan mode ..." << endl;
      F.generateFragments(TermExtension::COMPLEXITY_SCAN);
    } else {
      cout << "Generating fragments using center residue mode..." << endl;
      F.generateFragments(TermExtension::CEN_RES);
    }
    timer.stop();
    cout << "Took " << timer.getDuration() << "s to generate fragments" << endl;
    
    ////Extend fragments
    timer.start();
    cout << "Extending the fragments using the MANY_CONTACT mode..." << endl;
    string extfragDir_aa = extfragDir + aa + "_";
    F.extendFragmentsandWriteStructures(seedTERM::seedType::MANY_CONTACT, extfragDir_aa);
    timer.stop();
    cout << "Took " << timer.getDuration() << "s to extend fragments" << endl;
    cout << "Generated " << F.getExtendedFragmentNumber() << " extended fragments" << endl;
    fs << aa << "\t" << F.getExtendedFragmentNumber() << endl;
    
    F.resetFragments();
    
    ////Write seeds to file for visualization
    string line_file = extfragDir + aa + "_seeds.lcloud";
    fstream line_fs;
    MstUtils::openFile(line_fs,line_file,fstream::out);
    
    StructuresBinaryFile seeds(extfragDir_aa+"extendedfragments.bin");
    while (seeds.hasNext() == true) {
      Structure* seed = seeds.next();
      Chain* seed_C = seed->getChainByID("0");
      classifier.writeCaInfotoLineFile(seed_C, line_fs);
    }
  }
  
  // Write the fragment structures/information
  string fragDir = "./";
  cout << "Storing fragments in '" << fragDir << "'..." << endl;
  F.writeFragmentPDBs(fragDir);
  
  return 0;
}
