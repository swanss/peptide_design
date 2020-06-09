//
//  peptideWindowDistance.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 6/5/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include "mstoptions.h"
#include "mstrotlib.h"
#include "mstsystem.h"

#include "mstsystem_exts.h"
#include "utilities.h"
#include "benchmarkutilities.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Takes a list of peptide-protein structures, creates backbone windows from the peptides, records their distance from the protein");
  op.addOption("list", "The path to a list of pdb files. Each structure must contain at least a peptide and a protein, with a specifically formatted name (e.g. 5UUL_A_B)",true);
  op.addOption("config", "The path to a configuration file",true);
  op.setOptions(argc, argv);
  
  //define variables
  int flanking_residues = 2;
  mstreal cd_cut = 0.1;
  mstreal int_cut = 0.1;
  mstreal bb_dist = 3.5;
  
  string list_file = op.getString("list");
  vector<string> pdb_paths = MstUtils::fileToArray(list_file);
  
  configFile config(op.getString("config"));
  
  fstream info_file;
  MstUtils::openFile(info_file,"./peptide_window_distance.tsv");
  info_file << "name\tpeptide_chain_id\twindow\tdistance" << endl;
  
  for (string pdb_path : pdb_paths) {
    string pdb_name = MstSystemExtension::fileName(pdb_path);
    
    vector<string> split = MstUtils::split(pdb_name);
    string peptide_ID = split[2];
    
    Structure complex(pdb_path);
    Chain* p_chain = complex.getChainByID(peptide_ID);
    vector<Residue*> protein_residue;
    for (Residue* R: p_chain->getResidues()) {
      protein_residue.push_back(R);
    }

    seedStatistics seedStat(complex,peptide_ID);
    
    Structure peptide(*p_chain);
    Structure peptide_backbone = RotamerLibrary::getBackbone(peptide);
    
    ConFind C(config.getRL(),complex);
    
    vector<pair<Residue*, Residue*>> contacts = generalUtilities::getContactsWith(protein_residue, C, 0, cd_cut, int_cut, bb_dist, true);
    vector<Residue*> peptide_contacts = generalUtilities::getContactingResidues(contacts);
    
    for (Residue* R: peptide_contacts) {
      //get the residue in the peptide_backbone structure
      Residue* R_pep = peptide_backbone[0].findResidue(R->getName(),R->getNum());
      int i = R_pep->getResidueIndex();
      int n_i = max(0,i-flanking_residues);
      int c_i = min(i+flanking_residues,peptide_backbone.residueSize());
      
      //construct windows around each peptide residue that contacts the protein
      AtomPointerVector window;
      for (int id = n_i; id <= c_i; id++) {
        Residue* R_window = &peptide_backbone.getResidue(id);
        window.push_back(R_window);
      }
      
      //compute the distance of each window from the protein and write to the file
      CartesianPoint center = window.getGeometricCenter();
      mstreal distance = seedStat.point2NearestProteinAtom(center);
      
      info_file << pdb_name << "\t" << p_chain->getID() << "\t" << n_i << "\t" << distance << endl;
    }
  }
  info_file.close();
  return 0;
};
