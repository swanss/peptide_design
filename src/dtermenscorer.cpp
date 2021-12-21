#include "dtermenscorer.h"

dTERMenScorer::dTERMenScorer() {
  init();
}

dTERMenScorer::dTERMenScorer(const string& configFile) {
  init();
  readConfigFile(configFile);
}

void dTERMenScorer::init() {
  kT = 1.0;
  aaMapType = 1;
  cdCut = 0.01;
  intCut = 0.01; // set to a value over 1.0 to not count sidechain-backbone contacts via interference
  pmSelf = 1;
  pmPair = 1;
  selfResidualPC = selfCorrPC = 1.0;
  selfResidualMinN = 1000;
  selfResidualMaxN = 5000;
  selfCorrMaxCliqueSize = -1;
  selfCorrMinN = 200;
  selfCorrMaxN = 5000;
  pairMinN = 1000;
  pairMaxN = 5000;
  recordData = false;
  homCut = 0.6;
  setAminoAcidMap();
  setEnergyFunction("35");

  // set up FASST base options
  F.setOptions(fasstSearchOptions());
  F.setRedundancyProperty("sim");
  foptsBase = F.options();
}

void dTERMenScorer::setEnergyFunction(const string& ver) {
  efunVer = ver;
  if (efunVer.compare("35") == 0) {
    // this is the default
  } else if (efunVer.compare("35.f") == 0) {
    // do not use interference-based contacts
    intCut = 1.1;
  } else {
    MstUtils::error("unknown energy function version '" + efunVer + "'");
  }
}

void dTERMenScorer::readConfigFile(const string& configFile) {
  vector<string> lines = MstUtils::fileToArray(configFile);
  for (int i = 0; i < lines.size(); i++) {
    string line = MstUtils::trim(MstUtils::removeComment(lines[i], "#"));
    if (line.empty()) continue;
    vector<string> ents = MstUtils::trim(MstUtils::split(line, "="));
    if (ents.size() != 2) MstUtils::error("could not parse parameter line '" + lines[i] + "' from file " + configFile, "dTERMenScorer::dTERMen(const string&)");
    if (ents[0].compare("fasstdb") == 0) {
      fasstdbPath = ents[1];
    } else if (ents[0].compare("rotlib") == 0) {
      rotLibFile = ents[1];
    } else if (ents[0].compare("efun") == 0) {
      setEnergyFunction(ents[1]);
    } else if (ents[0].compare("selfCorrMaxCliqueSize") == 0) {
      selfCorrMaxCliqueSize = MstUtils::toInt(ents[1]);
    } else if (ents[0].compare("selfResidualLims") == 0) {
      vector<int> lims = MstUtils::splitToInt(ents[1]);
      if (lims.size() != 2) MstUtils::error("expected two integers in selfResidualLims field", "dTERMenScorer::readConfigFile");
      selfResidualMinN = lims[0];
      selfResidualMaxN = lims[1];
    } else if (ents[0].compare("pairLims") == 0) {
      vector<int> lims = MstUtils::splitToInt(ents[1]);
      if (lims.size() != 2) MstUtils::error("expected two integers in pairLims field", "dTERMenScorer::readConfigFile");
      pairMinN = lims[0];
      pairMaxN = lims[1];
    } else if (ents[0].compare("selfCorrLims") == 0) {
      vector<int> lims = MstUtils::splitToInt(ents[1]);
      if (lims.size() != 2) MstUtils::error("expected two integers in selfCorrLims field", "dTERMenScorer::readConfigFile");
      selfCorrMinN = lims[0];
      selfCorrMaxN = lims[1];
    } else if (ents[0].compare("homCut") == 0) {
      homCut = MstUtils::toReal(ents[1]);
    } else {
      MstUtils::error("unknown parameter name '" + ents[0] + "'", "dTERMenScorer::dTERMen(const string&)");
    }
  }

  if (fasstdbPath.empty()) MstUtils::error("FASST database not defined in configuration file " + configFile, "dTERMenScorer::dTERMen(const string&)");
  F.readDatabase(fasstdbPath, 2);
  if (!backPotFile.empty()) {
    readBackgroundPotentials(backPotFile);
  } else {
    buildBackgroundPotentials();
  }
  if (rotLibFile.empty()) { MstUtils::error("dTERMen configuration file does not specify a rotamer library, '" + configFile + "'", "dTERMenScorer::readConfigFile"); }
  else {
    RL.readRotamerLibrary(rotLibFile);
  }
}

vector<pair<Residue*, Residue*>> dTERMenScorer::getContactsWith(const vector<Residue*>& source, ConFind& C, int type, bool verbose) {
  set<Residue*> sourceSet = MstUtils::contents(source);
  vector<pair<Residue*, Residue*>> conts;
  set<pair<Residue*, Residue*>> contsSet;
  pair<Residue*, Residue*> c;
  contactList contList;
  if (verbose) {
    cout << "identifying contacts with:";
    for (int i = 0; i < source.size(); i++) cout << " " << *(source[i]);
    cout << endl;
  }

  // get all contacts involving the source residues
  for (int cType = 0; cType < 2; cType++) {
    if (cType == 0) contList = C.getContacts(source, cdCut);
    else {
      if (intCut > 1.0) continue;
      contList = C.getInterfering(source, intCut);
    }

    // go through each and insert into list, in the right order, if it qualifies
    contList.sortByDegree();
    for (int i = 0; i < contList.size(); i++) {
      bool isInA = (sourceSet.find(contList.residueA(i)) != sourceSet.end());
      bool isInB = (sourceSet.find(contList.residueB(i)) != sourceSet.end());
      if (((type == 0) && (isInA == isInB)) || ((type == 1) && !(isInA && isInB)) || ((type == 2) && !(isInA || isInB))) continue;

      if (!isInA) c = pair<Residue*, Residue*>(contList.residueB(i), contList.residueA(i));
      else c = pair<Residue*, Residue*>(contList.residueA(i), contList.residueB(i));

      if (verbose) cout << "\t" << (cType ? "interference " : "contact-degree ") << "contact with " << *(c.second) << " (from " << *(c.first) << "); " << contList.degree(i) << endl;
      if (contsSet.find(c) == contsSet.end()) {
        contsSet.insert(c);
        conts.push_back(c);
      }
    }
  }
  if (verbose) cout << "in the end, found " << conts.size() << " contacts" << endl;

  return conts;
}

void dTERMenScorer::setAminoAcidMap() {
  /* Perfectly corresponding to standard residues. */
  map<string, string> standard = {{"HSD", "HIS"}, {"HSE", "HIS"}, {"HSC", "HIS"}, {"HSP", "HIS"}};

  /* Almost perfectly corresponding to standard residues:
   * MSE -- selenomethyonine; SEC -- selenocysteine */
  map<string, string> almostStandard = {{"MSE", "MET"}, {"SEC", "CYS"}};

  /* A little less perfectly corresponding pairings, but can be acceptable (depends):
   * HIP -- ND1-phosphohistidine; SEP -- phosphoserine; TPO -- phosphothreonine;
   * PTR -- o-phosphotyrosine. */
  map<string, string> lessStandard = {{"HIP", "HIS"}, {"SEP", "SER"}, {"TPO", "THR"}, {"PTR", "TYR"}};

  for (res_t aai = 0; aai <= SeqTools::maxIndex(); aai++) {
    string aa = SeqTools::idxToTriple(aai);
    if (aai <= 19) {
      aaMap[aai] = aai; // natural amino acids go as they are, of course
    } else {
      switch (aaMapType) {
        case 1:
          // map frequent modifications to their closest standard amino acids
          // discard amino acids not explicitly mentioned in the above sets
          if (standard.find(aa) != standard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(standard[aa]);
          } else if (almostStandard.find(aa) != almostStandard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(almostStandard[aa]);
          } else if (lessStandard.find(aa) != lessStandard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(lessStandard[aa]);
          }
          break;
        case 2:
          // map only the obvious modifications to their closest standard amino acid,
          // but mostly preserve the modifications (good for designing with modifications)
          // discard amino acids not explicitly mentioned in the above sets
          if (standard.find(aa) != standard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(standard[aa]);
          } else if (almostStandard.find(aa) != almostStandard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(almostStandard[aa]);
          } else if (lessStandard.find(aa) != lessStandard.end()) {
            aaMap[aai] = aai;
          }
          break;
        case 3:
          // do no chemical mapping, interpret everything explicitly, BUT still
          // discard amino acids not explicitly mentioned in the above sets
          if (standard.find(aa) != standard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(standard[aa]);
          } else if (almostStandard.find(aa) != almostStandard.end()) {
            aaMap[aai] = aai;
          } else if (lessStandard.find(aa) != lessStandard.end()) {
            aaMap[aai] = aai;
          }
          break;
        case 4:
          // do no chemical mapping, interpret everything explicitly
          if (standard.find(aa) != standard.end()) {
            aaMap[aai] = SeqTools::aaToIdx(standard[aa]);
          } else {
            aaMap[aai] = aai;
          }
          break;
        default:
          MstUtils::error("unrecognized amino-acid mapping type " + MstUtils::toString(aaMapType));
      }
    }
  }
  for (res_t aai = 0; aai <= SeqTools::maxIndex(); aai++) {
    if (aaMap.find(aai) == aaMap.end()) continue; // not an allowed amino acid
    int i = aaMap[aai];
    if (aaIdx.find(i) == aaIdx.end()) { int idx = aaIdx.size(); aaIdx[i] = idx; }
  }
  globAlph.resize(aaIdx.size());
  for (auto it = aaIdx.begin(); it != aaIdx.end(); ++it) {
    globAlph[it->second] = it->first;
  }
}

void dTERMenScorer::printAminoAcidMap() {
  cout << "---- amino-acid map ----" << endl;
  cout << "universal alphabet has " << globAlph.size() << " amino acids:\n\t";
  for (int i = 0; i < globAlph.size(); i++) cout << SeqTools::idxToTriple(globAlph[i]) << " (" << aaIdx[globAlph[i]] << "); ";
  cout << endl << "allowed mappings are:" << endl;
  for (auto it = aaMap.begin(); it != aaMap.end(); ++it) {
    cout << SeqTools::idxToTriple(it->first) << " -> " << SeqTools::idxToTriple(it->second) << endl;
  }
}

bool dTERMenScorer::isInGlobalAlphabet(const string& aa) const {
  return (aaToIndex(aa) >= 0);
}

bool dTERMenScorer::isInGlobalAlphabet(res_t aa) const {
  return (aaToIndex(aa) >= 0);
}

int dTERMenScorer::aaToIndex(const string& aa) const {
  return aaToIndex(SeqTools::aaToIdx(aa));
}

int dTERMenScorer::aaToIndex(res_t aa) const {
  if (aaIdx.find(aa) != aaIdx.end()) return aaIdx.at(aa);
  return -1;
}

bool dTERMenScorer::aaIndexKnown(int aaIdx) const {
  return aaIdx >= 0;
}

string dTERMenScorer::indexToResName(int idx) const {
  return SeqTools::idxToTriple(globAlph[idx]);
}

res_t dTERMenScorer::indexToAA(int idx) const {
  return globAlph[idx];
}

void dTERMenScorer::buildBackgroundPotentials() {
  // extract all necessary residue properties
  vector<string> propNames = {"phi", "psi", "omega", "env"};
  map<string, vector<mstreal> > propVals;
  vector<int> aa;
  for (int i = 0; i < propNames.size(); i++) MstUtils::assert(F.isResiduePropertyDefined(propNames[i]), "property " + propNames[i] + " is not defined in the FASST database", "dTERMenScorer::buildBackgroundPotentials()");

  for (int ti = 0; ti < F.numTargets(); ti++) {
    Sequence S = F.getTargetSequence(ti);
    int N = S.length();

    // compute multiplicity of each residue
    vector<mstreal> mult(N, 1); // multiplicity of each residue in the structure
    map<int, vector<FASST::resAddress>> simsToTarget = F.getResidueRelationships(ti, "sim");
    for (auto it = simsToTarget.begin(); it != simsToTarget.end(); ++it) {
      mult[it->first] += (it->second).size();
    }

    // store all properties
    for (int ri = 0; ri < N; ri++) {
      string aaName = S.getResidue(ri, true);
      if (!isInGlobalAlphabet(aaName)) continue;
      aa.push_back(aaToIndex(aaName));
      for (int i = 0; i < propNames.size(); i++) {
        propVals[propNames[i]].push_back(F.getResidueProperty(ti, propNames[i], ri));
      }
      propVals["mult"].push_back(mult[ri]);
    }
  }

  // for each position in the database, accumulate total statistical potential,
  // for  all possible amino acids, from all known pseudo-energy types
  vector<vector<mstreal> > back(aa.size(), vector<mstreal> (globalAlphabetSize(), 0.0));
  bkPot = buildZeroDimPotential(aa, back);
  // cout << "Background frequency potential:\n"; printZeroDimPotential(bkPot);
  ppPot = buildTwoDimPotential(binData(propVals["phi"], propVals["psi"], {-180, 180, 36}, {-180, 180, 36}, propVals["mult"], true), aa, 10.0, back, true);
  // cout << "Phi/psi potential:\n"; printTwoDimPotential(ppPot);
  omPot = buildOneDimPotential(binData(propVals["omega"], 2, {-180, 180, 1000, 1}, propVals["mult"], true), aa, 10.0, back, true);
  // cout << "Omega potential:\n"; printOneDimPotential(omPot);
  envPot = buildOneDimPotential(binData(propVals["env"], 1, {0, 1, 75}, propVals["mult"], false), aa, 10.0, back, true);
  // cout << "Env potential:\n"; printOneDimPotential(envPot);
// TODO: could also do an end potential: 1, 2, 3 treated specially and N-2, N-1
// N also.
}

int dTERMenScorer::findBin(const vector<mstreal>& binEdges, mstreal x) {
  if ((x < binEdges.front()) || (x > binEdges.back())) return -1;
  // do a binary search
  int nb = binEdges.size() - 1;
  int left = 0, right = nb - 1, k;
  while (true) {
    k = (left + right)/2;
    if (x < binEdges[k]) {
      right = k - 1;
    } else if ((x > binEdges[k+1]) || ((k < nb - 1) && (x == binEdges[k+1]))) { // the last bin includes the opper limit
      left = k + 1;
    } else {
      break;
    }
  }
  return k;
}

mstreal dTERMenScorer::lookupZeroDimPotential(const zeroDimPotType& P, int aa) {
  return P.aaEnergies[aa];
}

mstreal dTERMenScorer::lookupOneDimPotential(const oneDimPotType& P, mstreal x, int aa) {
  int k = findBin(P.binEdges, x);
  if (k < 0) return 0;
  return P.aaEnergies[k][aa];
}

mstreal dTERMenScorer::lookupTwoDimPotential(const twoDimPotType& P, mstreal x, mstreal y, int aa) {
  int kx = findBin(P.xBinEdges, x);
  if (kx < 0) return 0;
  int ky = findBin(P.yBinEdges, y);
  if (ky < 0) return 0;
  return P.aaEnergies[kx][ky][aa];
}

dTERMenScorer::twoDimHist dTERMenScorer::binData(const vector<mstreal>& X, const vector<mstreal>& Y, const vector<mstreal>& xBinSpec, const vector<mstreal>& yBinSpec, const vector<mstreal>& M, bool isAngle) {
  dTERMenScorer::oneDimHist xH = binData(X, 1, xBinSpec, M, isAngle);
  dTERMenScorer::oneDimHist yH = binData(Y, 1, yBinSpec, M, isAngle);
  dTERMenScorer::twoDimHist xyH;
  xyH.xBinEdges = xH.binEdges;
  xyH.yBinEdges = yH.binEdges;
  xyH.weights = xH.weights; // the weights vector should be the same for X- and Y- histograms, so can copy from either
  xyH.bins.resize(xH.bins.size(), vector<vector<int> >(yH.bins.size()));

  // find intersections between X and Y bins
  map<int, int> yPointsToBins;
  for (int i = 0; i < yH.bins.size(); i++) {
    for (int k = 0; k < yH.bins[i].size(); k++) {
      yPointsToBins[yH.bins[i][k]] = i;
    }
  }
  for (int i = 0; i < xH.bins.size(); i++) {
    for (int k = 0; k < xH.bins[i].size(); k++) {
      int idx = xH.bins[i][k];
      if (yPointsToBins.find(idx) != yPointsToBins.end()) {
        int j = yPointsToBins[idx];
        xyH.bins[i][j].push_back(idx);
// cout << "[" << X[idx] << ", " << Y[idx] << "] -> " << "{" << xyH.xBinEdges[i] << ":" << xyH.xBinEdges[i+1] << ", " << xyH.yBinEdges[j] << ":" << xyH.yBinEdges[j+1] << "}" << endl;
      }
    }
  }

  return xyH;
}

dTERMenScorer::oneDimHist dTERMenScorer::binData(const vector<mstreal>& X, int binSpecType, const vector<mstreal>& binSpec, const vector<mstreal>& M, bool isAngle) {
  // limit to data points within range and transform dihedral angles
  mstreal minVal = binSpec[0];
  mstreal maxVal = binSpec[1];
  vector<mstreal> x = X, m = M;
  vector<int> origIndices = MstUtils::range(0, (int) X.size());
  vector<int> exclude;
  for (int i = 0; i < x.size(); i++) {
    if ((x[i] < minVal) || (x[i] > maxVal)) {
      exclude.push_back(i);
      continue;
    }
    if (isAngle) {
      mstreal an = CartesianGeometry::angleDiff(x[i], 0);
      // represent +/- pi as -pi; then, the bin boundary definition always works
      if (an >= 180) an = -180;
      x[i] = an;
    }
  }
  if (m.empty()) m.resize(x.size(), 1.0);
  if (!exclude.empty()) {
    vector<mstreal> cleanedData(x.size() - exclude.size());
    vector<mstreal> cleanedMult(cleanedData.size());
    vector<int> cleanedIndices(cleanedData.size());
    int k = 0, j = 0;
    for (int i = 0; i < x.size(); i++) {
      if (i == exclude[k]) { k++; continue; }
      cleanedData[j] = x[i];
      cleanedMult[j] = m[i];
      cleanedIndices[j] = origIndices[i];
      j++;
    }
    x = cleanedData;
    origIndices = cleanedIndices;
    m = cleanedMult;
  }

  // actually do binning
  vector<int> sortedInds = MstUtils::sortIndices(x);
  oneDimHist H;
  if (binSpecType == 1) {
    /* Uniform binning */
    if (binSpec.size() != 3) MstUtils::error("expected three values for bin specification type 1", "dTERMenScorer::binData");
    int nb = (int) binSpec[2];
    if ((nb <= 0) || (minVal >= maxVal)) MstUtils::error("inconsistent bin specification, type 1", "dTERMenScorer::binData");
    mstreal bw = (maxVal - minVal)/nb;
    H.binEdges.resize(nb + 1, 0);
    for (int i = 0; i < nb; i++) H.binEdges[i] = minVal + bw*i;
    H.binEdges[nb] = maxVal;
  } else if (binSpecType == 2) {
    /* Non-uniform binning with some minimal number of elements per bin and a
     * minimal bin width. */
    if (binSpec.size() != 4) MstUtils::error("expected four values for bin specification type 2", "dTERMenScorer::binData");
    int minNumPerBin = (int) binSpec[2];
    mstreal minBinWidth = binSpec[3];
    if ((minNumPerBin <= 0) || (minVal >= maxVal) || (minBinWidth < 0)) MstUtils::error("inconsistent bin specification, type 2", "dTERMenScorer::binData");
    if (minNumPerBin > x.size()) MstUtils::error("requested min number of elements per bin exceeds the total number of elements", "dTERMenScorer::binData");
    if (minBinWidth > maxVal - minVal) MstUtils::error("requested min bin width exceeds specified range", "dTERMenScorer::binData");
    /* To account for multiplicities when counting bin sizes, create a vector of
     * cumulative data "mass", whereby massI[k] stores the mass for the first k
     * smallest elements (i.e., for indices in range sortedInds[0], sortedInds[1],
     * ..., sortedInds[k]). massE[k] is the same, but it excludes the mass of
     * the final point sortedInds[k] (having both arrays simplifies some counting
     * later). The mass of each data point is 1 over its multiplicity (if it is
     * defined). Thus, if multiplicity is 1 (i.e., the point is "unique"), then
     * its mass is 1. But if multiplicity is 2, meaning there are two roughly
     * equivalent positions in the dataset, then each will together count as 1
     * (each counting as 1/2). If multiplicities are not given, then all counts
     * are set to 1, which means data counts are interpreted explicitly. */
    vector<mstreal> massI(x.size()), massE(x.size());
    if (m.size() != x.size()) MstUtils::error("number of points and multiplicities specified is not consistent!", "dTERMenScorer::binData");
    massI[0] = (1/m[sortedInds[0]]); massE[0] = 0;
    for (int i = 1; i < m.size(); i++) {
      massI[i] = massI[i-1] + 1/m[sortedInds[i]];
      massE[i] = massE[i-1] + 1/m[sortedInds[i-1]];
    }

    vector<mstreal> binEdgesFromLeft, binEdgesFromRight;
    binEdgesFromLeft.push_back(minVal);
    binEdgesFromRight.push_back(maxVal);
    // leftInd and rightInd will always store the first index (from left or right,
    // respectively) that maps into the next bin. Thus, by setting rightIndex to
    // the last index, the right edge of the last bin will be inclusive
    int leftInd = 0, rightInd = x.size() - 1;
    while (true) {
      int numRemaining = massI[rightInd] - massE[leftInd];
      mstreal remWidth = x[sortedInds[rightInd]] - x[sortedInds[leftInd]];
      if ((numRemaining >= 2*minNumPerBin) && (remWidth >= 2*minBinWidth)) { // if enough points left for at least two bins
        // add a bin on the left
        for (int k = leftInd; k < sortedInds.size(); k++) {
          if ((massE[k] - massI[leftInd] >= minNumPerBin) && (x[sortedInds[k]] - binEdgesFromLeft.back() > minBinWidth)) {
            binEdgesFromLeft.push_back(x[sortedInds[k]]);
            leftInd = k; // the left-most element in the next bin (left-to-right)
            break;
          }
        }

        // add a bin on the right
        for (int k = rightInd; k >= 0; k--) {
          if ((massI[rightInd] - massI[k] >= minNumPerBin) && (binEdgesFromRight.back() - x[sortedInds[k]] > minBinWidth)) {
            binEdgesFromRight.push_back(x[sortedInds[k]]);
            rightInd = k - 1; // the right-most element in the next bin (right-to-left)
            break;
          }
        }

        // if there are enough points left for no more than two bins, just split
        // the rest between the two bins
        if ((numRemaining < 3*minNumPerBin) || (remWidth < 3*minBinWidth)) {
          if (rightInd < leftInd) {
            /* Even though in this iteration there were initially enough points to fit three bins
            * (by width and number), it is possible that to add a left and right bin (respecting
            * minimal width and number of points) can move the terminal indices past each other.
            * Example, sorted values between leftInd and rightInd are: [0.5, 0.51, 0.52, 0.53,
            * 0.54, 0.55, 0.6, 0.7] minimum bin width is 0.1, and minimum number of points per bin
            * is 3. Then, to satisfy all conditions going from the left we need to do [0.5 - 0.6),
            * but to satisfy all conditions going from the right would require [0.54 - 0.7). The
            * two segments overlap. If this happens, we will make just one bin between leftInd and
            * rightInd, since it is not possible to make two while respecting all rules. */
            binEdgesFromLeft.pop_back();
            binEdgesFromRight.pop_back();
            H.binEdges = binEdgesFromLeft;
            H.binEdges.insert(H.binEdges.end(), binEdgesFromRight.rbegin(), binEdgesFromRight.rend());
          } else {
            int k = (leftInd + rightInd)/2;
            H.binEdges = binEdgesFromLeft;
            H.binEdges.back() = x[sortedInds[k]];
            H.binEdges.insert(H.binEdges.end(), binEdgesFromRight.rbegin() + 1, binEdgesFromRight.rend());
          }
          break;
        }
      } else if ((numRemaining >= minNumPerBin) && (remWidth >= minBinWidth)) { // if only enough points left for one bin
        // add the last bin to bridge the gap
        H.binEdges = binEdgesFromLeft;
        H.binEdges.insert(H.binEdges.end(), binEdgesFromRight.rbegin(), binEdgesFromRight.rend());
        break;
      } else {
        MstUtils::error("this should not have happened! Ran out of points/width for the middle bin.", "dTERMenScorer::binData");
      }
    }
  }

  // classify points into bins
  H.bins.resize(H.binEdges.size() - 1);
  H.binMasses.resize(H.binEdges.size() - 1);
  H.weights.resize(X.size(), 0.0);
  for (int i = 0; i < X.size(); i++) { // copy all weights (even for excluded points)
    H.weights[i] = M.empty() ? 1.0 : 1.0/M[i];
  }
  int bi = 0; mstreal binMass = 0;
  for (int i = 0; i < sortedInds.size(); i++) {
    // the last bin includes the right-most end of the range
    if ((x[sortedInds[i]] >= H.binEdges[bi+1]) && ((bi < H.bins.size() - 1) || (x[sortedInds[i]] > H.binEdges[bi+1]))) {
      H.binMasses[bi] = binMass;
      binMass = 0;
      bi++;
    }
    H.bins[bi].push_back(origIndices[sortedInds[i]]);
    binMass += 1/m[sortedInds[i]];
  }

  return H;
}

dTERMenScorer::zeroDimPotType dTERMenScorer::buildZeroDimPotential(const vector<int>& AA, vector<vector<mstreal> >& backPot) {
  // get the potential
  zeroDimPotType pot;
  int naa = globalAlphabetSize();
  map<int, mstreal> aaCounts;
  for (int i = 0; i < naa; i++) aaCounts[i] = 0.0;
  for (int i = 0; i < AA.size(); i++) aaCounts[AA[i]] += 1;
  pot.aaEnergies = vector<mstreal>(naa, 1.0/0.0);
  for (int i = 0; i < naa; i++) {
    pot.aaEnergies[i] = -kT*log(aaCounts[i]/AA.size());
  }

  // compute the background energy
  for (int i = 0; i < AA.size(); i++) {
    for (int j = 0; j < naa; j++) {
      backPot[i][j] = pot.aaEnergies[j];
    }
  }

  return pot;
}

dTERMenScorer::oneDimPotType dTERMenScorer::buildOneDimPotential(const oneDimHist& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot) {
  oneDimPotType pot;
  pot.binEdges = H.binEdges;
  int nb = pot.binEdges.size() - 1;
  int naa = globalAlphabetSize();
  pot.aaEnergies.resize(nb, vector<mstreal> (naa, 0.0));

  // overall amino-acid frequencies
  CartesianPoint fAA(naa);
  for (int bi = 0; bi < nb; bi++) { // iterating over bin indices ignores any points excluded from binning
    const vector<int>& binInds = H.bins[bi];
    for (int i = 0; i < binInds.size(); i++) {
      int k = binInds[i];
      fAA[AA[k]] += H.weights[k];
    }
  }
  fAA /= fAA.sum();

  // add up the weighted counts of each amino acid in each bin
  for (int bi = 0; bi < nb; bi++) {
    vector<mstreal> expCounts = vector<mstreal> (naa, 0.0); // expectation of each amino acid in this bin
    const vector<int>& binInds = H.bins[bi];
    for (int i = 0; i < binInds.size(); i++) {
      int k = binInds[i];
      mstreal w = H.weights[k]; // position weight
      int aa = AA[k];
      pot.aaEnergies[bi][aa] += w;

      // compute the expectation of each amino acid in this position
      if (!backPot.empty()) {
        mstreal Z = 0; // partition function for this position (over all amino acids)
        for (int aai = 0; aai < naa; aai++) Z += exp(-backPot[k][aai]/kT);
        for (int aai = 0; aai < naa; aai++) expCounts[aai] += w*exp(-backPot[k][aai]/kT)/Z;
      } else {
        // if no prior potential given, assume uniform prior
        for (int aai = 0; aai < naa; aai++) expCounts[aai] += w/naa;
      }
    }

    // compute the potential using frequency-proportional pseudocounts
    for (int aa = 0; aa < naa; aa++) {
      pot.aaEnergies[bi][aa] = -kT*log((pot.aaEnergies[bi][aa] + pc*fAA[aa]*naa)/(expCounts[aa] + pc*fAA[aa]*naa));
    }
  }

  // update the prior potential
  if (!backPot.empty() && updateBackPot) {
    for (int bi = 0; bi < nb; bi++) {
      const vector<int>& binInds = H.bins[bi];
      for (int i = 0; i < binInds.size(); i++) {
        int k = binInds[i];
        for (int aai = 0; aai < naa; aai++) {
          backPot[k][aai] += pot.aaEnergies[bi][aai];
        }
      }
    }
  }

  return pot;
}


dTERMenScorer::twoDimPotType dTERMenScorer::buildTwoDimPotential(const twoDimHist& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot) {
  twoDimPotType pot;
  pot.xBinEdges = H.xBinEdges;
  pot.yBinEdges = H.yBinEdges;
  int xnb = H.xBinEdges.size() - 1;
  int ynb = H.yBinEdges.size() - 1;
  int naa = globalAlphabetSize();
  pot.aaEnergies.resize(xnb, vector<vector<mstreal> >(ynb, vector<mstreal>(naa, 0.0)));

  // overall amino-acid frequencies
  CartesianPoint fAA(naa);
  for (int bi = 0; bi < xnb; bi++) { // iterating over bin indices ignores any points excluded from binning
    for (int bj = 0; bj < ynb; bj++) {
      const vector<int>& binInds = H.bins[bi][bj];
      for (int i = 0; i < binInds.size(); i++) {
        int k = binInds[i];
        fAA[AA[k]] += H.weights[k];
      }
    }
  }
  fAA /= fAA.sum();

  // add up the weighted counts of each amino acid in each bin
  for (int bi = 0; bi < xnb; bi++) { // iterating over bin indices ignores any points excluded from binning
    for (int bj = 0; bj < ynb; bj++) {
      vector<mstreal> expCounts = vector<mstreal> (naa, 0.0); // expectation of each amino acid in this bin
      const vector<int>& binInds = H.bins[bi][bj];
      for (int i = 0; i < binInds.size(); i++) {
        int k = binInds[i];
        mstreal w = H.weights[k]; // position weight
        int aa = AA[k];
        pot.aaEnergies[bi][bj][aa] += w;

        // compute the expectation of each amino acid in this position
        if (!backPot.empty()) {
          mstreal Z = 0; // partition function for this position (over all amino acids)
          for (int aai = 0; aai < naa; aai++) Z += exp(-backPot[k][aai]/kT);
          for (int aai = 0; aai < naa; aai++) expCounts[aai] += w*exp(-backPot[k][aai]/kT)/Z;
        } else {
          // if no prior potential given, assume uniform prior
          for (int aai = 0; aai < naa; aai++) expCounts[aai] += w/naa;
        }
      }

      // compute the potential using frequency-proportional pseudocounts
      for (int aa = 0; aa < naa; aa++) {
        pot.aaEnergies[bi][bj][aa] = -kT*log((pot.aaEnergies[bi][bj][aa] + pc*fAA[aa]*naa)/(expCounts[aa] + pc*fAA[aa]*naa));
      }
    }
  }

  // update the prior potential
  if (!backPot.empty() && updateBackPot) {
    for (int bi = 0; bi < xnb; bi++) { // iterating over bin indices ignores any points excluded from binning
      for (int bj = 0; bj < ynb; bj++) {
        const vector<int>& binInds = H.bins[bi][bj];
        for (int i = 0; i < binInds.size(); i++) {
          int k = binInds[i];
          for (int aai = 0; aai < naa; aai++) {
            backPot[k][aai] += pot.aaEnergies[bi][bj][aai];
          }
        }
      }
    }
  }

  return pot;
}

dTERMenScorer::oneDimPotType dTERMenScorer::buildOneDimPotential(const oneDimHist& H, const vector<int>& AA, mstreal pc) {
  vector<vector<mstreal> > backPot;
  return buildOneDimPotential(H, AA, pc, backPot, false);
}

void dTERMenScorer::printZeroDimPotential(const zeroDimPotType& P) {
  int naa = P.aaEnergies.size();
  for (int aa = 0; aa < naa; aa++) {
    cout << indexToResName(aa) << "\t" << P.aaEnergies[aa] << endl;
  }
}

void dTERMenScorer::printOneDimPotential(const oneDimPotType& P) {
  int nb = P.aaEnergies.size();
  if (nb == 0) return;
  int naa = globalAlphabetSize();
  for (int aa = 0; aa < naa; aa++) cout << indexToResName(aa) << " ";
  cout << endl;
  for (int bi = 0; bi < nb; bi++) {
    printf("%f %f", P.binEdges[bi], P.binEdges[bi+1]);
    for (int aa = 0; aa < naa; aa++) {
      printf(" %5.3f", P.aaEnergies[bi][aa]);
    }
    printf("\n");
  }
}

void dTERMenScorer::printTwoDimPotential(const twoDimPotType& P) {
  int xnb = P.xBinEdges.size() - 1;
  int ynb = P.yBinEdges.size() - 1;
  int naa = globalAlphabetSize();
  for (int aa = 0; aa < naa; aa++) {
    cout << "---------------> " << indexToResName(aa) << ":" << endl;
    for (int i = 0; i < xnb; i++) {
      for (int j = 0; j < ynb; j++) {
        printf("%5.3f ", P.aaEnergies[i][j][aa]);
      }
      printf("\n");
    }
  }
}

void dTERMenScorer::readBackgroundPotentials(const string& file) {

}

void dTERMenScorer::writeRecordedData(const string& file) {
  // find any duplicate TERMs
  map<set<int>, int> byResidueSet;
  for (int i = 0; i < data.size(); i++) {
    set<int> resIdxSet = MstUtils::contents(data[i].getResidueIndices());
    if (byResidueSet.find(resIdxSet) == byResidueSet.end()) {
      byResidueSet[resIdxSet] = i;
    } else {
      // if a TERM with the same residues was seen before, keep the one with more matches
      int pi = byResidueSet[resIdxSet];
      if (data[i].numMatches() > data[pi].numMatches()) byResidueSet[resIdxSet] = i;
    }
  }
  vector<int> unique = MstUtils::values(byResidueSet);
  sort(unique.begin(), unique.end());

  // write data about the target
  fstream ofs; MstUtils::openFile(ofs, file, ios::out);
  ofs << targetOrigSeq.toString() << endl;
  ofs << MstUtils::vecToString(variableResidues) << endl;
  for (int i = 0; i < targetOrigSeq.length(); i++) {
    ofs << targetResidueProperties["phi"][i] << " ";
    ofs << targetResidueProperties["psi"][i] << " ";
    ofs << targetResidueProperties["omega"][i] << " ";
    ofs << targetResidueProperties["env"][i] << endl;
  }

  // write data about unique TERMs
  for (int i : unique) {
    ofs << "* TERM " << i << endl;
    ofs << MstUtils::vecToString(data[i].getResidueIndices()) << endl;
    for (int j = 0; j < data[i].numMatches(); j++) {
      const fasstSolution& m = data[i].getMatch(j);
      ofs << F.getMatchSequence(m).toString() << " " << m.getRMSD() << " ";
      ofs << MstUtils::vecToString(F.getResidueProperties(m, "phi")) << " ";
      ofs << MstUtils::vecToString(F.getResidueProperties(m, "psi")) << " ";
      ofs << MstUtils::vecToString(F.getResidueProperties(m, "omega")) << " ";
      ofs << MstUtils::vecToString(F.getResidueProperties(m, "env")) << endl;
    }
  }
  ofs.close();
}

mstreal dTERMenScorer::selfEnergy(Residue* R, const string& aa) {
  return (selfEnergies(R))[dTERMenScorer::aaToIndex(aa.empty() ? R->getName() : aa)];
}

vector<mstreal> dTERMenScorer::selfEnergies(Residue* R, bool verbose) {
  if (R->getStructure() == NULL) MstUtils::error("cannot operate on a disembodied residue!", "dTERMenScorer::selfEnergies(Residue*, bool)");
  ConFind C(&RL, *(R->getStructure()));
  return selfEnergies(R, C, verbose);
}

vector<mstreal> dTERMenScorer::selfEnergies(Residue* R, ConFind& C, bool verbose, bool seedIndicator) {
  auto rmsdCutSelfRes = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.0, 20); };
  auto rmsdCutSelfCor = [](const vector<int>& fragResIdx, const Structure& S) { return RMSDCalculator::rmsdCutoff(fragResIdx, S, 1.1, 15); };
  if (R->getStructure() == NULL) MstUtils::error("cannot operate on a disembodied residue!", "dTERMenScorer::selfEnergies(Residue*, ConFind&, bool)");
  Structure& S = *(R->getStructure());

  // -- simple environment components
  if (verbose) cout << "\tdTERMenScorer::selfEnergies -> trivial statistical components..." << endl;
  int naa = globalAlphabetSize();
  CartesianPoint selfE(naa, 0.0);
  for (int aai = 0; aai < naa; aai++) {
    selfE[aai] = backEner(aai) + bbOmegaEner(R->getOmega(), aai) + bbPhiPsiEner(R->getPhi(), R->getPsi(), aai) + envEner(C.getFreedom(R), aai);
  }
  if (verbose) printSelfComponent(selfE, "\t");

  // -- self residual
  if (verbose) cout << "\tdTERMenScorer::selfEnergies -> self residual..." << endl;
  termData sT({R}, pmSelf);
  F.setOptions(foptsBase);
  F.setQuery(sT.getTERM());
  F.setRMSDCutoff(rmsdCutSelfRes(sT.getResidueIndices(), S));
  F.setMinNumMatches(selfResidualMinN);
  F.setMaxNumMatches(selfResidualMaxN);
  sT.setMatches(F.search(), homCut, &F);
  CartesianPoint selfResidual = singleBodyStatEnergy(sT.getMatches(), sT.getCentralResidueIndices()[0], selfResidualPC);
  if (recordData) data.push_back(sT);
  if (verbose) printSelfComponent(selfResidual, "\t");
  selfE += selfResidual;

  if (!seedIndicator) return selfE;

  // -- self correction
  if ((selfCorrMaxCliqueSize >= 0) && (selfCorrMaxCliqueSize < 2)) return selfE; // if max clique size is less than 2, then there is effectively no self residual
  if (verbose) cout << "\tdTERMenScorer::selfEnergies -> self correction..." << endl;

  // -- get contacts
  vector<pair<Residue*, Residue*>> conts = getContactsWith({R}, C, 0, verbose);
  vector<Residue*> contResidues(conts.size(), NULL);
  for (int i = 0; i < conts.size(); i++) contResidues[i] = conts[i].second;

  // consider each contacting residue with the central one and see if there are
  // enough matches. If so, create a clique to be grown later. If not, still
  // create a clique (with number of matches), which will not be grown.
  vector<termData> finalCliques;
  map<Residue*, termData> cliquesToGrow;
  for (int i = 0; i < contResidues.size(); i++) {
    if (verbose) cout << "\t\tdTERMenScorer::selfEnergies -> seed clique with contact " << *(contResidues[i]) << "..." << endl;
    termData c({R, contResidues[i]}, pmSelf);
    F.setOptions(foptsBase);
    F.setQuery(c.getTERM());
    F.setRMSDCutoff(rmsdCutSelfCor(c.getResidueIndices(), S));
    F.setMinNumMatches(selfCorrMinN);
    F.setMaxNumMatches(selfCorrMaxN);
    c.setMatches(F.search(), homCut, &F);
    if ((c.numMatches() < selfCorrMinN) || (c.getMatch(selfCorrMinN - 1).getRMSD()) > F.getRMSDCutoff()) { finalCliques.push_back(c); }
    else { cliquesToGrow[contResidues[i]] = c; }
  }

  // for those contacts with sufficient matches, try to combine into larger cliques
  while (!cliquesToGrow.empty()) {
    // sort remaining contacting residues by decreasing number of matches of the
    // clique comprising the residue and the central residue
    vector<Residue*> remConts = MstUtils::keys(cliquesToGrow);
    sort(remConts.begin(), remConts.end(), [&cliquesToGrow](Residue* i, Residue* j) { return cliquesToGrow[i].numMatches() > cliquesToGrow[j].numMatches(); });

    // pick the one with highest number of contacts, and try to grow the clique
    termData parentClique = cliquesToGrow[remConts[0]];
    remConts.erase(remConts.begin());
    termData grownClique = parentClique;
    if (verbose) cout << "\t\tdTERMenScorer::selfEnergies -> will try to grow [" << parentClique.toString() << "]..." << endl;
    while (!remConts.empty()) {
      if ((selfCorrMaxCliqueSize >= 0) && (parentClique.numCentralResidues() >= selfCorrMaxCliqueSize)) break;
      // try to add every remaining contact
      for (int j = 0; j < remConts.size(); j++) {
        if (verbose) cout << "\t\t\tdTERMenScorer::selfEnergies -> trying to add " << *(remConts[j]) << "..." << endl;
        termData newClique = parentClique;
        newClique.addCentralResidue(remConts[j], pmSelf);
        F.setOptions(foptsBase);
        F.setQuery(newClique.getTERM());
        F.setRMSDCutoff(rmsdCutSelfCor(newClique.getResidueIndices(), S));
        F.setMaxNumMatches(selfCorrMaxN);
        newClique.setMatches(F.search(), homCut, &F);
        if ((j == 0) || (newClique.numMatches() > grownClique.numMatches())) {
          if (verbose) cout << "\t\t\t\tdTERMenScorer::selfEnergies -> new best" << endl;
          grownClique = newClique;
        }
      }
      // EITHER a sufficient number of matches OR not too many fewer than for the
      // parent clique is the logic we had implemented in original dTERMen
      if ((grownClique.numMatches() >= selfCorrMinN) || (grownClique.numMatches() > 0.8*parentClique.numMatches())) {
        remConts = MstUtils::setdiff(remConts, {grownClique.getCentralResidues().back()});
        parentClique = grownClique;
        if (verbose) cout << "\t\t\tdTERMenScorer::selfEnergies -> chose to add " << *(grownClique.getCentralResidues().back()) << "..." << endl;
      } else {
        if (verbose) cout << "\t\t\tdTERMenScorer::selfEnergies -> nothing more worked, sticking with parent clique..." << endl;
        grownClique = parentClique;
        break;
      }
    }

    // add the grown clique to the list of final cliques, and remove used up
    // residues from the list of cliques to grow
    if (verbose) cout << "\t\tdTERMenScorer::selfEnergies -> adding final clique [" << grownClique.toString() << "]..." << endl;
    finalCliques.push_back(grownClique);
    vector<Residue*> deleted = MstUtils::setdiff(MstUtils::keys(cliquesToGrow), remConts);
    if (verbose) cout << "\t\tdTERMenScorer::selfEnergies -> neighboring residues taken care of this cycle:";
    for (int i = 0; i < deleted.size(); i++) {
      if (verbose) cout << " " << *(deleted[i]);
      cliquesToGrow.erase(deleted[i]);
    }
    if (verbose) cout << endl;
  }

  // finally, actually compute the self correction for each final clique
  if (verbose) cout << "\tdTERMenScorer::selfEnergies -> final cliques:" << endl;
  for (int i = 0; i < finalCliques.size(); i++) {
    if (verbose) cout << "\t\t" << finalCliques[i].toString() << endl;
    if (recordData) data.push_back(finalCliques[i]);
    CartesianPoint cliqueDelta = singleBodyStatEnergy(finalCliques[i].getMatches(), finalCliques[i].getCentralResidueIndices()[0], selfCorrPC);
    if (verbose) printSelfComponent(cliqueDelta, "\t\t\t");
    selfE += cliqueDelta;
  }

  return selfE;
}

void dTERMenScorer::printSelfComponent(const CartesianPoint& ener, const string& prefix) {
  cout << prefix;
  for (int i = 0; i < globAlph.size(); i++) printf("%8s", indexToResName(i).c_str());
  cout << endl << prefix;
  for (int i = 0; i < ener.size(); i++) printf("%8.3f", ener[i]);
  cout << endl;
}

CartesianPoint dTERMenScorer::singleBodyStatEnergy(fasstSolutionSet& matches, int cInd, int pc) {
  CartesianPoint selfE(globalAlphabetSize(), 0.0);
  CartesianPoint Ne = dTERMenScorer::singleBodyExpectations(matches, cInd);
  CartesianPoint No = dTERMenScorer::singleBodyObservations(matches, cInd);
  for (int aai = 0; aai < selfE.size(); aai++) {
    selfE[aai] = -kT*log((No[aai] + pc)/(Ne[aai] + pc));
  }

  return selfE;
}

CartesianPoint dTERMenScorer::singleBodyObservations(fasstSolutionSet& matches, int cInd) {
  CartesianPoint No(globalAlphabetSize(), 0.0);
  for (int i = 0; i < matches.size(); i++) {
    int aaIdx = dTERMenScorer::aaToIndex(F.getMatchSequence(matches[i])[cInd]);
    if (!aaIndexKnown(aaIdx)) continue; // this match has some amino acid that is not known in the current alphabet
    No[aaIdx] += 1.0;
  }
  return No;
}

CartesianPoint dTERMenScorer::twoBodyObservations(fasstSolutionSet& matches, int cIndI, int cIndJ) {
  CartesianPoint No(globalAlphabetSize()*globalAlphabetSize(), 0.0);
  for (int i = 0; i < matches.size(); i++) {
    int aaiIdx = dTERMenScorer::aaToIndex(F.getMatchSequence(matches[i])[cIndI]);
    int aajIdx = dTERMenScorer::aaToIndex(F.getMatchSequence(matches[i])[cIndJ]);
    if (!aaIndexKnown(aaiIdx) || !aaIndexKnown(aajIdx)) continue;
    No[dTERMenScorer::pairToIdx(aaiIdx, aajIdx)] += 1.0;
  }
  return No;
}

CartesianPoint dTERMenScorer::singleBodyExpectations(fasstSolutionSet& matches, int cInd, vector<CartesianPoint>* breakDown) {
  CartesianPoint Ne(globalAlphabetSize(), 0.0);
  if (breakDown != NULL) { breakDown->clear(); breakDown->resize(matches.size()); }
  for (int i = 0; i < matches.size(); i++) {
    int aaIdx = dTERMenScorer::aaToIndex(F.getMatchSequence(matches[i])[cInd]);
    if (!aaIndexKnown(aaIdx)) continue; // this match has some amino acid that is not known in the current alphabet
    CartesianPoint inMatchExp = backExpectation(matches[i], cInd);
    if (breakDown != NULL) (*breakDown)[i] = inMatchExp;
    Ne += inMatchExp;
  }
  return Ne;
}

vector<CartesianPoint> dTERMenScorer::singleBodyExpectationsMatchedMarginals(fasstSolutionSet& matches, int cInd) {
  // compute initial expected amino-acid distributions in each match
  vector<CartesianPoint> mlogP;
  CartesianPoint Ne(globalAlphabetSize(), 0.0);
  vector<int> validMatches;
  for (int i = 0; i < matches.size(); i++) {
    int aaIdx = dTERMenScorer::aaToIndex(F.getMatchSequence(matches[i])[cInd]);
    if (!aaIndexKnown(aaIdx)) continue; // this match has some amino acid that is not known in the current alphabet
    validMatches.push_back(i);
    CartesianPoint inMatchExp = backExpectation(matches[i], cInd);
    Ne += inMatchExp;
    mlogP.push_back(inMatchExp);
    for (int j = 0; j < mlogP.back().size(); j++) mlogP.back()[j] = -log(mlogP.back()[j]);
  }
  CartesianPoint No = dTERMenScorer::singleBodyObservations(matches, cInd);

  // next, iterate to find optimal amino-acid bias energies to get marginals to agree
  CartesianPoint delE(globalAlphabetSize(), 0.0);
  vector<CartesianPoint> breakDown(matches.size());
  int Ncyc = 10;
  for (int c = 0; c < Ncyc; c++) {
    for (int aai = 0; aai < delE.size(); aai++) delE[aai] = delE[aai] + log(Ne[aai]/No[aai]);
    Ne.clear(); Ne.resize(globalAlphabetSize(), 0.0);
    for (int i = 0; i < mlogP.size(); i++) {
      CartesianPoint inMatchExp = mlogP[i] + delE;
      dTERMenScorer::enerToProb(inMatchExp);
      Ne += inMatchExp;
      if (c == Ncyc - 1) {
        breakDown[validMatches[i]] = inMatchExp;
      }
    }
  }

  return breakDown;
}

CartesianPoint dTERMenScorer::backExpectation(const fasstSolution& m, int cInd) {
  vector<mstreal> p(globalAlphabetSize(), 0);
  mstreal phi = F.getResidueProperties(m, "phi")[cInd];
  mstreal psi = F.getResidueProperties(m, "psi")[cInd];
  mstreal omg = F.getResidueProperties(m, "omega")[cInd];
  mstreal env = F.getResidueProperties(m, "env")[cInd];
  for (int aai = 0; aai < globalAlphabetSize(); aai++) {
    p[aai] = backEner(aai) + bbOmegaEner(omg, aai) + bbPhiPsiEner(phi, psi, aai) + envEner(env, aai);
  }
  dTERMenScorer::enerToProb(p);
  return p;
}

mstreal dTERMenScorer::enerToProb(vector<mstreal>& ener) {
  mstreal minEner = MstUtils::min(ener);
  mstreal Z = 0;
  for (int i = 0; i < ener.size(); i++) { ener[i] = exp(-(ener[i] - minEner)/kT); Z += ener[i]; }
  for (int i = 0; i < ener.size(); i++) ener[i] /= Z;
  return Z;
}

int dTERMenScorer::pairToIdx(int aai, int aaj) const {
  return aai*globalAlphabetSize() + aaj;
}

pair<int, int> dTERMenScorer::idxToPair(int idx) const {
  return pair<int, int>(idx / globalAlphabetSize(), idx % globalAlphabetSize());
}