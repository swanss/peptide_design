//
//  testConFind.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 11/11/19.
//

#include "mstsystem.h"
#include "mstoptions.h"
#include "mstcondeg.h"

#include "utilities.h"

int main(int argc, char *argv[]) {
  MstOptions op;
  op.setTitle("For testing the different types of contacts and the tolerateMissingBackbone option in ConFind.");
  op.addOption("p", "The structure which will be used to compute contact information",true);
  op.addOption("mode","The type of contacts that will be computed. Options: cont, int, bbInt.", true);
  op.addOption("config", "The path to the configuration file",true);
  op.addOption("sel","A selection string that specifies for which residues to look for interactions.");
  op.addOption("unk_sel","A selection string that specifies which residue should have names set to 'UNK'");
  op.setOptions(argc,argv);
  
  configFile config(op.getString("config"));
  string RL = config.getRL();

  Structure S(op.getString("p"));
  ConFind C(RL,S,true);
  
  if (op.isGiven("unk_sel")) {
    selector sel(S);
    vector<Residue*> selection = sel.selectRes(op.getString("unk_sel"));
    cout << "Selected " << selection.size() << " residues to set to 'UNK'" << endl;
    for (Residue* R : selection) R->setName("UNK");
  }
  
  vector<Residue*> selection;
  if (op.isGiven("sel")) {
    selector sel(S);
    selection = sel.selectRes(op.getString("sel"));
  }
  
  string contact_type = op.getString("mode");
  cout << "Contact type: " << contact_type << " selected" << endl;
  
  contactList list;
  mstreal threshold = .01;
  mstreal dist = 3.5;
  
  if (contact_type == "cont") {
    if (!op.isGiven("sel")) list = C.getContacts(S,threshold);
    else list = C.getContacts(selection,threshold);
  } else if (contact_type == "int") {
    if (!op.isGiven("sel")) list = C.getInterference(S,threshold);
    else list = C.getInterference(selection,threshold);
  } else if (contact_type == "bbInt") {
    if (!op.isGiven("sel"))list = C.getBBInteraction(S,dist);
    else list = C.getBBInteraction(selection,dist);
  } else {
    MstUtils::error("Contact type is not recognized!");
  }
  
  fstream out;
  string output_path = "./"+ contact_type + ".txt";
  MstUtils::openFile(out, output_path, fstream::out);
  
  list.sortByDegree();
  cout << "Now writing contacts..." << endl;
  for (int i = 0; i < list.size(); i++) {
    out << list.residueA(i)->getChainID() << list.residueA(i)->getNum() << "\t" << list.residueB(i)->getChainID() << list.residueB(i)->getNum() << "\t" << list.degree(i) << endl;
  }
  out.close();
  
  //get the unique destination contacts
  vector<Residue*> dest = list.destResidues();
  vector<Residue*> dest_unique;
  for (Residue* R : dest) if (find(dest_unique.begin(),dest_unique.end(),R) == dest_unique.end()) dest_unique.push_back(R);
  
  output_path = "./dest_conts_selstring_"+ contact_type + ".txt";
  MstUtils::openFile(out, output_path, fstream::out);
  out << generalUtilities::selectionStringFromRes(dest_unique) << endl;
}
