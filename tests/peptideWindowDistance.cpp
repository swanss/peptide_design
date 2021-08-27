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
  mstreal cd_cut = 0.01;
  mstreal int_cut = 0.01;
  mstreal bb_dist = 3.5;
  
  string list_file = op.getString("list");
  vector<string> pdb_paths = MstUtils::fileToArray(list_file);
  cout << "list had " << pdb_paths.size() << " structures" << endl;
  
  configFile config(op.getString("config"));
  
  fstream info_file;
  MstUtils::openFile(info_file,"peptide_window_distance.tsv",fstream::out);
  info_file << "name\tpeptide_chain_id\twindow\tdistance" << endl;
  
  for (string pdb_path : pdb_paths) {
    string pdb_name = MstSystemExtension::fileName(pdb_path);
    cout << "pdb file name: " << pdb_name << endl;
    
    vector<string> split = MstUtils::split(pdb_name,"_");
    string peptide_ID = split[2];
    cout << "peptide chain name: " << peptide_ID << endl;
    
    //load the complex and construct a ConFind object
    Structure complex(pdb_path);
    Structure cleaned_complex;
    RotamerLibrary::extractProtein(cleaned_complex,complex);
    ConFind C(config.getRL(),cleaned_complex);
    
    //get the protein residues and find all peptide residues that contact these
    Chain* pep_chain = cleaned_complex.getChainByID(peptide_ID);
    vector<Residue*> protein_residue;
    for (int i = 0; i < cleaned_complex.chainSize(); i++) {
      Chain* C = &cleaned_complex.getChain(i);
      if (C != pep_chain) {
        vector<Residue*> chain_res = C->getResidues();
        protein_residue.insert(protein_residue.end(),chain_res.begin(),chain_res.end());
      }
    }
    
    vector<pair<Residue*, Residue*>> contacts = generalUtilities::getContactsWith(protein_residue, C, 0, cd_cut, int_cut, bb_dist, true);
    vector<Residue*> peptide_contacts = generalUtilities::getContactingResidues(contacts);
    
    //get the peptide backbone (for constructing windows that are comparable to the seeds)
    Structure peptide(*pep_chain);
    Structure peptide_backbone(RotamerLibrary::getBackbone(peptide));
    cout << "got backbone of peptide" << endl;
    
    seedStatistics seedStat(cleaned_complex,peptide_ID);
    for (Residue* R: peptide_contacts) {
      //get the residue in the peptide_backbone structure
      Residue* R_pep = peptide_backbone[0].findResidue(R->getName(),R->getNum());
      int i = R_pep->getResidueIndex();
      int n_i = max(0,i-flanking_residues);
      int c_i = min(i+flanking_residues,peptide_backbone.residueSize()-1); //rip, you never grow out of off-by-one errors, do you?
      
      //construct windows around each peptide residue that contacts the protein
      AtomPointerVector window;
      for (int id = n_i; id <= c_i; id++) {
        Residue* R_window = &peptide_backbone.getResidue(id);
        window.push_back(R_window);
      }
      
      //compute the distance of each window from the protein and write to the file
      CartesianPoint center = window.getGeometricCenter();
      mstreal distance = seedStat.point2NearestProteinAtom(center);
      
      info_file << pdb_name << "\t" << pep_chain->getID() << "\t" << n_i << "\t" << distance << endl;
    }
  }
  info_file.close();
  return 0;
};
