//
//  testFuseFragments.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 8/13/19.
//


#include "msttypes.h"
#include "mstoptions.h"
#include "mstfuser.h"

using namespace std;
using namespace MST;

// Various snippets of code were copied from TERMify on 19/08/13

mstreal getRadius(const Structure& S) {
  selector sel(S);
  AtomPointerVector atoms = sel.select("name N or name CA or name C or name O");
  mstreal rad = 0;
  for (int i = 0; i < atoms.size(); i++) {
    for (int j = i+1; j < atoms.size(); j++) {
      mstreal d = atoms[i]->distance(atoms[j]);
      if (d > rad) rad = d;
    }
  }
  return rad;
}

AtomPointerVector getBackbone(const Structure& S, vector<int> residues) {
  AtomPointerVector atoms;
  vector<string> bba = {"N", "CA", "C", "O"};
  for (int i = 0; i < residues.size(); i++) {
    Residue& res = S.getResidue(residues[i]);
    for (int j = 0; j < bba.size(); j++) {
      atoms.push_back(res.findAtom(bba[j]));
    }
  }
  return atoms;
}

int main (int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("Fuse a set of overlapping fragments. Note: the residues of each fragment must already be renumbered according to the fusion topology.");
  op.addOption("l","File specifying each fragment pdb file, one fragment per line.",true);
  op.addOption("path","Path to the directory containing the fragments, e.g., /path/to/dir", true);
  op.addOption("base","The base name of the new PDB",true);
  op.addOption("f","If known, a structure that contains a fixed region of the future structure can be imported.",true);
  op.addOption("sel","The chain ID of the peptide chain. This is not added to the fusion topology. The rest of structure is included and treated as fixed.",true);
  op.setOptions(argc,argv);
  
  //Import the list and all PDBs in the list
  vector<string> all_fragment_paths;
  all_fragment_paths = MstUtils::fileToArray(op.getString("l"));
  
  vector<Structure> all_fragments;
  for (int i = 0; i < all_fragment_paths.size(); i++) {
    Structure frag(op.getString("path") + "/" + all_fragment_paths[i]);
    all_fragments.push_back(frag);
  }
  
  //Import the fixed structure
  Structure fixed_whole(op.getString("f"));
  
  //select the residues that should be included as fixed to the fusion topology
  selector sel(fixed_whole);
  string p_cid = op.getString("sel");
  vector<Residue*> fixed_residues = sel.selectRes("not chain "+p_cid);
  vector<Residue*> unfixed_residues = sel.selectRes("chain "+p_cid);
  Structure fixed_subset(fixed_residues);
  
  cout << "Selection: " << op.getString("sel") << " selected " << MstUtils::toString(fixed_subset.residueSize()) << " residues" << endl;
  
  /* It is possible that at some of the levels of complexity either 1) no fragment could be defined
   containing some region or 2) a fragment was defined, but there was no match in the database within
   the lRMSD threshold. In this case there would be positions in the fusion topology that lack atoms.
   As an easy fix, I will quickly read through the fragment residue numbers to identify the longest
   possible peptide designable by the given fragments. Following this, I will adjust the numbers/
   the overall topology length accordingly.
   
   NOTE: This assumes that there is a path between each of the peptide fragments.
   NOTE: this assumes that the peptide follows the protein in the PDB of the complex */
  
  // find the range of the peptide
  int last_idx = fixed_whole.residueSize() - 1;
  int first_idx = fixed_subset.residueSize(); //technically + 1 - 1
  
  // find the highest and lowest residue number within this range
  int min = 10000;
  int max = 0;
  for (int frag = 0; frag < all_fragments.size(); frag++) {
    for (int res = 0; res < all_fragments[frag].residueSize(); res++) {
      Residue* R = &(all_fragments[frag].getResidue(res));
      int r_num = R->getNum();
      if ((r_num <= last_idx) && (r_num >= first_idx)) {
        if (r_num < min) {
          min = r_num;
        }
        if (r_num > max) {
          max = r_num;
        }
      }
    }
  }
  
  // renumber each fragment according to the lowest residue number
  for (int frag = 0; frag < all_fragments.size(); frag++) {
    for (int res = 0; res < all_fragments[frag].residueSize(); res++) {
      Residue* R = &(all_fragments[frag].getResidue(res));
      int r_num = R->getNum();
      
      //only renumber if in the range of the peptide
      if (r_num >= first_idx) R->setNum(r_num - min + first_idx);
    }
  }
  
  
  ////Build the fusionTopology object
  // make the fusion topology object using the adjusted highest number residue
  max = first_idx + max - min + 1; //last_idx after adjustment + 1 is total res size
  fusionTopology top(max);
  
  //Build vector of fixed residue indices
  vector<int> fixed;
  map<Residue*, int> indices = MstUtils::indexMap(fixed_whole.getResidues());
  fixed.resize(fixed_residues.size());
  for (int i = 0; i < fixed_residues.size(); i++) fixed[i] = indices[fixed_residues[i]];
  
  //Build a vector of the unfixed (peptide) residue indices
  vector<int> unfixed; unfixed.resize(unfixed_residues.size());
  for (int i = 0; i < unfixed_residues.size(); i++) unfixed[i] = indices[unfixed_residues[i]];
  
  //Set the fixed positions
  top.addFixedPositions(fixed);
  
  //Add the subset of the fixed structure
  top.addFragment(fixed_subset,fixed,-10.0); //by definition, fixed should also be able to specify the residueIDs (note
  
  //Add the fragments that contain new peptide residues
  for (int i = 0; i < all_fragments.size(); i++) {
    top.addFragment(all_fragments[i]);
  }
  
  // fuser options
  fusionParams opts; opts.setVerbose(false);
  opts.setMinimizerType(fusionParams::gradDescent);
  opts.setRepFC(1);
  opts.setCompFC(0.1);
  mstreal compactnessRadius = getRadius(fixed_whole);
  opts.setCompRad(compactnessRadius);
  cout << "will be trying to combine structure to a radius of " << compactnessRadius << endl;
  fusionOutput propScore;
  
  //Here I assume that all fragments must be included in the topology
  Structure fused = Fuser::fuse(top, propScore, opts);
  
  RMSDCalculator rc;
  
  //Realign to the fixed section for visualization
  AtomPointerVector before = getBackbone(fixed_whole, fixed);
  AtomPointerVector after = getBackbone(fused, fixed);
  rc.align(after, before, fused);
  
  //Calculate the RMSD over the design to the native peptide
  //get the native peptide chain
  AtomPointerVector native_peptide = getBackbone(fixed_whole, unfixed);
  
  //get the fused peptide chain
  AtomPointerVector fused_peptide = getBackbone(fused, unfixed);
  
  mstreal rmsd = rc.rmsd(native_peptide,fused_peptide);
  string rmsd_name = MstUtils::toString(rmsd);
  
  cout << "the RMSD between the native and fused peptide is " << rmsd_name << endl;
  
  //Save structure
  fused.writePDB(op.getString("base")+"_"+rmsd_name+".pdb");
};

