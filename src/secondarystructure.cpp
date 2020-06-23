//
//  secondarystructure.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 2/16/20.
//

#include "secondarystructure.h"

dihedralProbabilities::dihedralProbabilities(SecStructType type) {
  // 18 x 18 vector of vectors copied from rdmap.cpp in the STRIDE project folder
  // Rows correspond to Phi, Columns correspond to Psi
  /*
   -180  {{X, X, X},
    |       . . .
Phi 0     {X, X, X},
    |       . . .
   -180   {X, X, X}}
   
       -180 - 0 - 180
             Psi
   */
  if (type == SecStructType::HELIX) {
    probabilities = {{0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0009014423, 0.0041898815, 0.0085105160, 0.0133839026, 0.0245425366, 0.0407802090, 0.0464176536, 0.0330946408, 0.0134803243, 0.0024038462, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0007370283, 0.0077203326, 0.0269849468, 0.0492307022, 0.0621860325, 0.0747849122, 0.0919913873, 0.0918549150, 0.0617070347, 0.0241584498, 0.0041428790, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0041416897, 0.0287234355, 0.0835687742, 0.1384727061, 0.1562444866, 0.1470608264, 0.1360232681, 0.1159155145, 0.0742164999, 0.0290896539, 0.0050673936, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0009375000, 0.0156580955, 0.0757770315, 0.1856354773, 0.2785892785, 0.2880102694, 0.2332847565, 0.1741978228, 0.1281246394, 0.0793832615, 0.0320557840, 0.0058840578, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0048893229, 0.0437000208, 0.1617751122, 0.3399706185, 0.4626395404, 0.4418565035, 0.3235570788, 0.2100441158, 0.1358627081, 0.0776144490, 0.0297011137, 0.0052390974, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0136979166, 0.0917820632, 0.2773087323, 0.5047551394, 0.6214492917, 0.5485223532, 0.3655386865, 0.2054343373, 0.1121114418, 0.0548815951, 0.0178668182, 0.0025975490, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0246484373, 0.1396044195, 0.3594934344, 0.5710113049, 0.6337110400, 0.5133636594, 0.3054708838, 0.1402616948, 0.0584463216, 0.0228670351, 0.0058531328, 0.0005151099, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0265885405, 0.1365883052, 0.3163702190, 0.4545661211, 0.4628692269, 0.3425511420, 0.1761947423, 0.0607788190, 0.0158569515, 0.0042061093, 0.0008107311, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0152018229, 0.0738445148, 0.1630392224, 0.2269553691, 0.2237145752, 0.1528334022, 0.0652616471, 0.0150429625, 0.0014589608, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0035156249, 0.0165251363, 0.0379281938, 0.0584417619, 0.0619409233, 0.0404052660, 0.0136552500, 0.0016678370, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0011718750, 0.0046875002, 0.0070312503, 0.0046875002, 0.0011718750, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0006944445, 0.0036063762, 0.0080820229, 0.0101532144, 0.0076146079, 0.0032324446, 0.0006009616, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000}};
  }
  else if (type == SecStructType::SHEET) {
    probabilities = {{0.2769023776, 0.1408346891, 0.0464910716, 0.0073784725, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0047086575, 0.0218229108, 0.0569166169, 0.1254088134, 0.2340224832, 0.3511219919, 0.4355685711, 0.4584180117, 0.4007356465},
                    {0.4067636132, 0.2329865396, 0.0927943364, 0.0237838365, 0.0055147060, 0.0013786765, 0.0000000000, 0.0000000000, 0.0000000000, 0.0088186050, 0.0420726910, 0.1043856740, 0.2086037844, 0.3677131534, 0.5367187858, 0.6412357688, 0.6458424330, 0.5580080152},
                    {0.4286311865, 0.2678007782, 0.1282834113, 0.0529448465, 0.0220588241, 0.0055147060, 0.0000000000, 0.0000000000, 0.0000000000, 0.0086062262, 0.0445192643, 0.1197573245, 0.2487278134, 0.4369854629, 0.6241853237, 0.7160459757, 0.6829043031, 0.5716546178},
                    {0.3639202416, 0.2397334576, 0.1305907220, 0.0683420748, 0.0330882370, 0.0082720593, 0.0000000000, 0.0000000000, 0.0000000000, 0.0053559211, 0.0328565054, 0.1048930883, 0.2402425259, 0.4295993447, 0.6026929021, 0.6669865251, 0.6039550304, 0.4841639400},
                    {0.2637948096, 0.1723874062, 0.0920098722, 0.0464194641, 0.0220588241, 0.0055147060, 0.0000000000, 0.0000000000, 0.0000000000, 0.0030202419, 0.0224239044, 0.0804052502, 0.1923188865, 0.3456886411, 0.4811576009, 0.5223571062, 0.4586051404, 0.3565762639},
                    {0.1628032923, 0.0930610597, 0.0400134660, 0.0143100554, 0.0055147060, 0.0013786765, 0.0000000000, 0.0000000000, 0.0000000000, 0.0015453297, 0.0132468110, 0.0489843786, 0.1174781919, 0.2150468081, 0.3082944453, 0.3439011276, 0.3080393970, 0.2371628135},
                    {0.0825822726, 0.0338854715, 0.0092895878, 0.0012122844, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0003863324, 0.0046614520, 0.0186656341, 0.0477515720, 0.0961741805, 0.1546680480, 0.1961039603, 0.1944279373, 0.1469529718},
                    {0.0326442868, 0.0073916214, 0.0008854167, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0006347656, 0.0031504754, 0.0104655549, 0.0272454955, 0.0570511036, 0.0941907763, 0.1088592261, 0.0785619915},
                    {0.0090501504, 0.0007651417, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0007207961, 0.0035958111, 0.0131648667, 0.0318824202, 0.0425693691, 0.0292618107},
                    {0.0013020834, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0013020834, 0.0052083335, 0.0078125000, 0.0052083335},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
                    {0.0210939310, 0.0078523019, 0.0013020834, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0014204546, 0.0071634995, 0.0169352461, 0.0272206441, 0.0357281528, 0.0395361669, 0.0343801714},
                    {0.1146211401, 0.0503530800, 0.0130920913, 0.0015190972, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0010016026, 0.0046167620, 0.0157516468, 0.0453012958, 0.0937970504, 0.1454590708, 0.1861637682, 0.2019522935, 0.1764564067}};
  }
}

mstreal dihedralProbabilities::getResidueProb(Residue* R) {
  bool strict = false;
  mstreal phi = R->getPhi(strict);
  mstreal psi = R->getPsi(strict);
  //just return 0 if the dihedral can't be computed
  if (R->isBadDihedral(phi) || R->isBadDihedral(psi)) return 0.0;
  int phi_bin = angle2Bin(phi);
  int psi_bin = angle2Bin(psi);
  return probabilities[phi_bin][psi_bin];
}

int dihedralProbabilities::angle2Bin(mstreal angle) {
  return (angle + 180) / 20;
}

string secondaryStructureClassifier::classifyResidue(Residue* R) {
  mstreal helix_prob = helix.getResidueProb(R);
  mstreal sheet_prob = sheet.getResidueProb(R);
  
  string letter = (helix_prob > sheet_prob ? "H" : "E");
  dihedralProbabilities& prob = (helix_prob > sheet_prob ? helix : sheet);
  
  int same_class = 0;
  int max_flank = (2*flanking_res)+1;
  int actual_flank = max_flank;
  for (int offset = -flanking_res; offset <= flanking_res; offset++) {
    Residue* R_flanking = R->iPlusDelta(offset);
    if (R_flanking == NULL) actual_flank--;
    else if (prob.getResidueProb(R_flanking) > 0) {
      same_class++;
    }
  }
  float fraction_same_class = float(same_class)/float(actual_flank);
    if (fraction_same_class >= threshold) {
      return letter;
    } else return "O";
}


string secondaryStructureClassifier::classifyStruct(const Structure& S) {
  string classification = "";
  
  vector<Residue*> all_res = S.getResidues();
  
  for (Residue* R : all_res) {
    if (R->getAtom(0).isHetero()) continue;
    classification += classifyResidue(R);
  }
  
  return classification;
}

string secondaryStructureClassifier::classifyResInStruct(Structure* S, vector<int> all_res_idx) {
  string classification;
  for (int res_idx : all_res_idx) classification += classifyResidue(&S->getResidue(res_idx));
  return classification;
}

tuple<mstreal,mstreal,mstreal> secondaryStructureClassifier::getSecStructFractions(const Structure& Str) {
  int H = 0; int S = 0; int C = 0; //for counting observations of each class
  mstreal H_frac; mstreal S_frac; mstreal C_frac;
  
  vector<Residue*> residues = Str.getResidues();
  for (Residue* R : residues) {
    string classification = classifyResidue(R);
    if (classification == "H") {
      H += 1;
    } else if (classification == "E") {
      S += 1;
    } else if (classification == "O") {
      C += 1;
    }
  }
  H_frac = mstreal(H)/residues.size();
  S_frac = mstreal(S)/residues.size();
  C_frac = mstreal(C)/residues.size();
  return make_tuple(H_frac,S_frac,C_frac);
}

void secondaryStructureClassifier::writeResClassificationtoFile(Chain* C, fstream &out) {
  string structure_name = C->getStructure()->getName();
  
  out << structure_name;
  
  vector<Residue*> chain_res = C->getResidues();
  for (Residue* R : chain_res) {
    
    out << "," << classifyResidue(R);
  }
  out << endl;
}

void secondaryStructureClassifier::writeResClassificationtoTSV(Structure& S, fstream& out) {
  out << "chain_id\tresidue_num\tclassification" << endl;
  
  vector<Residue*> structure_res = S.getResidues();
  for (Residue* R : structure_res) {
    out << R->getChainID() << "\t";
    out << R->getNum() << "\t";
    out << classifyResidue(R);
    out << endl;
  }
}

void secondaryStructureClassifier::writeCentroidtoPointFile(string bin_path, int num_final_seeds, fstream& out, int bin_version) {
  StructuresBinaryFile bin_file(bin_path,true,bin_version);
  
  long num_seeds = bin_file.structureCount();
  mstreal skip_probability = 1 - min(mstreal(1),mstreal(num_final_seeds)/num_seeds);
  cout << "There are " << num_seeds << " seeds in the input binary file. Skip probability is " << skip_probability << endl;
  
  bin_file.reset(); int count = 0;
  while (bin_file.hasNext() == true) {
    mstreal sampled_value = MstUtils::randUnit();
    if (sampled_value <= skip_probability) {
      bin_file.skip();
      continue;
    }
    Structure* seed = bin_file.next();
    Chain* seed_C = seed->getChainByID("0");
    CartesianPoint centroid = AtomPointerVector(seed_C->getAtoms()).getGeometricCenter();
    out << centroid.getX() << "," << centroid.getY() << "," << centroid.getZ() << "," << 0 << endl;
    delete seed;
    count++;
  }
  out.close();
  cout << "In the end, wrote the coordinates of " << count << " seeds" << endl;
}

void secondaryStructureClassifier::writeCaInfotoPointFile(Chain* C, fstream &out) {
  vector<Residue*> chain_res = C->getResidues();
  for (Residue* R : chain_res) {
    Atom* Ca = R->findAtom("CA");
    out << Ca->getX() << "," << Ca->getY() << "," << Ca->getZ() << "," << classification2ColorID(classifyResidue(R)) << endl;
  }
}

void secondaryStructureClassifier::writeCaInfotoPointFile(string bin_path, fstream& out, int bin_version) {
  StructuresBinaryFile bin_file(bin_path,true,bin_version);
  bin_file.reset();
  while (bin_file.hasNext() == true) {
    Structure* seed = bin_file.next();
    Chain* seed_C = seed->getChainByID("0");
    writeCaInfotoPointFile(seed_C, out);
    delete seed;
  }
}

void secondaryStructureClassifier::writeCaInfotoLineFile(Chain* C, fstream& out) {
  vector<Residue*> chain_res = C->getResidues();
  for (Residue* R : chain_res) {
    Atom* Ca = R->findAtom("CA");
    string x_coord = MstUtils::toString(Ca->getX());
    string y_coord = MstUtils::toString(Ca->getY());
    string z_coord = MstUtils::toString(Ca->getZ());
    string classification = MstUtils::toString(classification2ColorID(classifyResidue(R)));
    string line = MstUtils::join(",",{x_coord,y_coord,z_coord,classification});
    if (R->nextResidue() != NULL) line += ",";
    out << line;
  }
  out << endl;
}

void secondaryStructureClassifier::writeCaInfotoLineFile(string bin_path, int num_final_seeds, fstream& out, int bin_version) {
  StructuresBinaryFile bin_file(bin_path,true,bin_version);
  
  long num_seeds = bin_file.structureCount();
  mstreal skip_probability = 1 - min(mstreal(1),mstreal(num_final_seeds)/num_seeds);
  cout << "There are " << num_seeds << " seeds in the input binary file. Skip probability is " << skip_probability << endl;
  
  bin_file.reset(); int count = 0;
  while (bin_file.hasNext() == true) {
    mstreal sampled_value = MstUtils::randUnit();
    if (sampled_value <= skip_probability) {
      bin_file.skip();
      continue;
    }
    Structure* seed = bin_file.next();
    Chain* seed_C = seed->getChainByID("0");
    writeCaInfotoLineFile(seed_C, out);
    delete seed;
    count++;
  }
  cout << "In the end, wrote the coordinates of " << count << " seeds" << endl;
  
}

void secondaryStructureClassifier::writeCaInfotoLineFile(string bin_path, vector<string> seed_names, fstream& out, int bin_version) {
  StructuresBinaryFile bin_file(bin_path,true,bin_version);
  bin_file.scanFilePositions();
  bin_file.reset();
  for (string seed_name : seed_names) {
    Structure* seed = bin_file.getStructureNamed(seed_name);
    Chain* seed_C = seed->getChainByID("0");
    writeCaInfotoLineFile(seed_C, out);
    delete seed;
  }
}
