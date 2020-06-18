#include "utilities.h"
#include "vdwRadii.h"

// bool vdwRadii::initConstants;
// double vdwRadii::maxSumRadii;
// double vdwRadii::getRadii;
// bool vdwRadii::contact;
// bool vdwRadii::independent;
// atomInteraction interactionType;
map<string, double> vdwRadii::radii;
bool vdwRadiiInitialized = vdwRadii::initConstants();


map<string, map<string, double>> resVDWRadii::radii;
bool resVDWRadiiInitialized = resVDWRadii::initConstants();

// add in residues
bool vdwRadii::initConstants() { radii = {{"N",  1.60}, {"NT", 1.60}, {"CA", 2.365}, {"C", 2.10}, {"CT", 2.10}, {"O", 1.60}, {"CB", 2.2350}}; return true; }
// bool vdwRadii::initConstants() { radii = {{"N", 1.85}, {"NT", 1.85}, {"CA", 2.275}, {"C", 2.0}, {"CT", 2.0}, {"O", 1.70}}; }


double vdwRadii::maxSumRadii() {
	double m = 0.0;
	for (auto it = radii.begin(); it != radii.end(); it++) {
		if(it->second > m) {
			m = it->second;
		}
	}
	return m * 2;
}

double vdwRadii::getRadii(const string& atomName) {
	return radii[atomName];
}

double vdwRadii::getRadii(const Atom& a) {
	return radii[a.getName()];
}

double vdwRadii::sumRadii(const Atom& a1, const Atom& a2) {
	return (radii[a1.getName()] + radii[a2.getName()]);
}

bool vdwRadii::clash(const Atom& a1, const Atom& a2, double lb) {
	return (a1.distance(a2) < vdwRadii::sumRadii(a1, a2) * lb);
}

bool vdwRadii::contact(const Atom& a1, const Atom& a2, double lb, double ub) {
	double dist = a1.distance(a2);
	double s = vdwRadii::sumRadii(a1, a2);
	return (dist >= s * lb && dist < s * ub);
}

bool vdwRadii::independent(const Atom& a1, const Atom& a2, double ub) {
	return (a1.distance(a2) >= vdwRadii::sumRadii(a1, a2) * ub);
}

atomInteraction vdwRadii::interactionType(const Atom& a1, const Atom& a2, double lb, double ub) {
	double dist = a1.distance(a2);
	double s = vdwRadii::sumRadii(a1, a2);
	if (dist < s * lb) {
		return CLASH;
	} else if (dist < s * ub) {
		return CONTACT;
	} else {
		return INDEPENDENT;
	}
}


// add in residues
bool resVDWRadii::initConstants() {
	radii["ALA"]["N"] = 1.6000;
	radii["ALA"]["H"] = 0.8000;
	radii["ALA"]["CA"] = 2.3650;
	radii["ALA"]["CB"] = 2.1650;
	radii["ALA"]["C"] = 2.1000;
	radii["ALA"]["O"] = 1.6000;
	radii["ARG"]["N"] = 1.6000;
	radii["ARG"]["H"] = 0.8000;
	radii["ARG"]["CA"] = 2.3650;
	radii["ARG"]["CB"] = 2.2350;
	radii["ARG"]["CG"] = 2.2350;
	radii["ARG"]["CD"] = 2.2350;
	radii["ARG"]["NE"] = 1.6000;
	radii["ARG"]["HE"] = 0.8000;
	radii["ARG"]["CZ"] = 2.1000;
	radii["ARG"]["NH1"] = 1.6000;
	radii["ARG"]["HH11"] = 0.6000;
	radii["ARG"]["HH12"] = 0.6000;
	radii["ARG"]["NH2"] = 1.6000;
	radii["ARG"]["HH21"] = 0.6000;
	radii["ARG"]["HH22"] = 0.6000;
	radii["ARG"]["C"] = 2.1000;
	radii["ARG"]["O"] = 1.6000;
	radii["ASN"]["N"] = 1.6000;
	radii["ASN"]["H"] = 0.8000;
	radii["ASN"]["CA"] = 2.3650;
	radii["ASN"]["CB"] = 2.2350;
	radii["ASN"]["CG"] = 2.1000;
	radii["ASN"]["OD1"] = 1.6000;
	radii["ASN"]["ND2"] = 1.6000;
	radii["ASN"]["HD21"] = 0.8000;
	radii["ASN"]["HD22"] = 0.8000;
	radii["ASN"]["C"] = 2.1000;
	radii["ASN"]["O"] = 1.6000;
	radii["ASP"]["N"] = 1.6000;
	radii["ASP"]["H"] = 0.8000;
	radii["ASP"]["CA"] = 2.3650;
	radii["ASP"]["CB"] = 2.2350;
	radii["ASP"]["CG"] = 2.1000;
	radii["ASP"]["OD1"] = 1.6000;
	radii["ASP"]["OD2"] = 1.6000;
	radii["ASP"]["C"] = 2.1000;
	radii["ASP"]["O"] = 1.6000;
	radii["AS1"]["N"] = 1.6000;
	radii["AS1"]["H"] = 0.8000;
	radii["AS1"]["CA"] = 2.3650;
	radii["AS1"]["CB"] = 2.2350;
	radii["AS1"]["CG"] = 2.1000;
	radii["AS1"]["OD2"] = 1.6000;
	radii["AS1"]["OD1"] = 1.6000;
	radii["AS1"]["HD1"] = 0.8000;
	radii["AS1"]["C"] = 2.1000;
	radii["AS1"]["O"] = 1.6000;
	radii["AS2"]["N"] = 1.6000;
	radii["AS2"]["H"] = 0.8000;
	radii["AS2"]["CA"] = 2.3650;
	radii["AS2"]["CB"] = 2.2350;
	radii["AS2"]["CG"] = 2.1000;
	radii["AS2"]["OD1"] = 1.6000;
	radii["AS2"]["OD2"] = 1.6000;
	radii["AS2"]["HD2"] = 0.8000;
	radii["AS2"]["C"] = 2.1000;
	radii["AS2"]["O"] = 1.6000;
	radii["CYS"]["N"] = 1.6000;
	radii["CYS"]["H"] = 0.8000;
	radii["CYS"]["CA"] = 2.3650;
	radii["CYS"]["CB"] = 2.2350;
	radii["CYS"]["SG"] = 1.8900;
	radii["CYS"]["C"] = 2.1000;
	radii["CYS"]["O"] = 1.6000;
	radii["GLN"]["N"] = 1.6000;
	radii["GLN"]["H"] = 0.8000;
	radii["GLN"]["CA"] = 2.3650;
	radii["GLN"]["CB"] = 2.2350;
	radii["GLN"]["CG"] = 2.2350;
	radii["GLN"]["CD"] = 2.1000;
	radii["GLN"]["OE1"] = 1.6000;
	radii["GLN"]["NE2"] = 1.6000;
	radii["GLN"]["HE21"] = 0.8000;
	radii["GLN"]["HE22"] = 0.8000;
	radii["GLN"]["C"] = 2.1000;
	radii["GLN"]["O"] = 1.6000;
	radii["GLU"]["N"] = 1.6000;
	radii["GLU"]["H"] = 0.8000;
	radii["GLU"]["CA"] = 2.3650;
	radii["GLU"]["CB"] = 2.2350;
	radii["GLU"]["CG"] = 2.2350;
	radii["GLU"]["CD"] = 2.1000;
	radii["GLU"]["OE1"] = 1.6000;
	radii["GLU"]["OE2"] = 1.6000;
	radii["GLU"]["C"] = 2.1000;
	radii["GLU"]["O"] = 1.6000;
	radii["GL1"]["N"] = 1.6000;
	radii["GL1"]["H"] = 0.8000;
	radii["GL1"]["CA"] = 2.3650;
	radii["GL1"]["CB"] = 2.2350;
	radii["GL1"]["CG"] = 2.2350;
	radii["GL1"]["CD"] = 2.1000;
	radii["GL1"]["OE2"] = 1.6000;
	radii["GL1"]["OE1"] = 1.6000;
	radii["GL1"]["HE1"] = 0.8000;
	radii["GL1"]["C"] = 2.1000;
	radii["GL1"]["O"] = 1.6000;
	radii["GL2"]["N"] = 1.6000;
	radii["GL2"]["H"] = 0.8000;
	radii["GL2"]["CA"] = 2.3650;
	radii["GL2"]["CB"] = 2.2350;
	radii["GL2"]["CG"] = 2.2350;
	radii["GL2"]["CD"] = 2.1000;
	radii["GL2"]["OE1"] = 1.6000;
	radii["GL2"]["OE2"] = 1.6000;
	radii["GL2"]["HE2"] = 0.8000;
	radii["GL2"]["C"] = 2.1000;
	radii["GL2"]["O"] = 1.6000;
	radii["GLY"]["N"] = 1.6000;
	radii["GLY"]["H"] = 0.8000;
	radii["GLY"]["CA"] = 2.2350;
	radii["GLY"]["C"] = 2.1000;
	radii["GLY"]["O"] = 1.6000;
	radii["HIS"]["N"] = 1.6000;
	radii["HIS"]["H"] = 0.8000;
	radii["HIS"]["CA"] = 2.3650;
	radii["HIS"]["CB"] = 2.2350;
	radii["HIS"]["CG"] = 2.1000;
	radii["HIS"]["ND1"] = 1.6000;
	radii["HIS"]["HD1"] = 0.8000;
	radii["HIS"]["CD2"] = 2.1000;
	radii["HIS"]["NE2"] = 1.6000;
	radii["HIS"]["CE1"] = 2.1000;
	radii["HIS"]["C"] = 2.1000;
	radii["HIS"]["O"] = 1.6000;
	radii["HSD"]["N"] = 1.6000;
	radii["HSD"]["H"] = 0.8000;
	radii["HSD"]["CA"] = 2.3650;
	radii["HSD"]["CB"] = 2.2350;
	radii["HSD"]["CG"] = 2.1000;
	radii["HSD"]["ND1"] = 1.6000;
	radii["HSD"]["CE1"] = 2.1000;
	radii["HSD"]["CD2"] = 2.1000;
	radii["HSD"]["NE2"] = 1.6000;
	radii["HSD"]["HE2"] = 0.8000;
	radii["HSD"]["C"] = 2.1000;
	radii["HSD"]["O"] = 1.6000;
	radii["HSC"]["N"] = 1.6000;
	radii["HSC"]["H"] = 0.8000;
	radii["HSC"]["CA"] = 2.3650;
	radii["HSC"]["CB"] = 2.2350;
	radii["HSC"]["CG"] = 2.1000;
	radii["HSC"]["CD2"] = 2.1000;
	radii["HSC"]["ND1"] = 1.6000;
	radii["HSC"]["HD1"] = 0.8000;
	radii["HSC"]["CE1"] = 2.1000;
	radii["HSC"]["NE2"] = 1.6000;
	radii["HSC"]["HE2"] = 0.8000;
	radii["HSC"]["C"] = 2.1000;
	radii["HSC"]["O"] = 1.6000;
	radii["ILE"]["N"] = 1.6000;
	radii["ILE"]["H"] = 0.8000;
	radii["ILE"]["CA"] = 2.3650;
	radii["ILE"]["CB"] = 2.3650;
	radii["ILE"]["CG2"] = 2.1650;
	radii["ILE"]["CG1"] = 2.2350;
	radii["ILE"]["CD1"] = 2.1650;
	radii["ILE"]["CD"] = 2.1650;
	radii["ILE"]["C"] = 2.1000;
	radii["ILE"]["O"] = 1.6000;
	radii["LEU"]["N"] = 1.6000;
	radii["LEU"]["H"] = 0.8000;
	radii["LEU"]["CA"] = 2.3650;
	radii["LEU"]["CB"] = 2.2350;
	radii["LEU"]["CG"] = 2.3650;
	radii["LEU"]["CD1"] = 2.1650;
	radii["LEU"]["CD2"] = 2.1650;
	radii["LEU"]["C"] = 2.1000;
	radii["LEU"]["O"] = 1.6000;
	radii["LYS"]["N"] = 1.6000;
	radii["LYS"]["H"] = 0.8000;
	radii["LYS"]["CA"] = 2.3650;
	radii["LYS"]["CB"] = 2.2350;
	radii["LYS"]["CG"] = 2.2350;
	radii["LYS"]["CD"] = 2.2350;
	radii["LYS"]["CE"] = 2.2350;
	radii["LYS"]["NZ"] = 1.6000;
	radii["LYS"]["HZ1"] = 0.6000;
	radii["LYS"]["HZ2"] = 0.6000;
	radii["LYS"]["HZ3"] = 0.6000;
	radii["LYS"]["C"] = 2.1000;
	radii["LYS"]["O"] = 1.6000;
	radii["MET"]["N"] = 1.6000;
	radii["MET"]["H"] = 0.8000;
	radii["MET"]["CA"] = 2.3650;
	radii["MET"]["CB"] = 2.2350;
	radii["MET"]["CG"] = 2.2350;
	radii["MET"]["SD"] = 1.8900;
	radii["MET"]["CE"] = 2.1650;
	radii["MET"]["C"] = 2.1000;
	radii["MET"]["O"] = 1.6000;
	radii["PEN"]["N"] = 1.6000;
	radii["PEN"]["H"] = 0.8000;
	radii["PEN"]["CA"] = 2.3650;
	radii["PEN"]["CB"] = 0.0000;
	radii["PEN"]["CG1"] = 2.1650;
	radii["PEN"]["CG2"] = 2.1650;
	radii["PEN"]["SG"] = 1.8900;
	radii["PEN"]["C"] = 2.1000;
	radii["PEN"]["O"] = 1.6000;
	radii["PHE"]["N"] = 1.6000;
	radii["PHE"]["H"] = 0.8000;
	radii["PHE"]["CA"] = 2.3650;
	radii["PHE"]["CB"] = 2.2350;
	radii["PHE"]["CG"] = 2.1000;
	radii["PHE"]["CD1"] = 2.1000;
	radii["PHE"]["CD2"] = 2.1000;
	radii["PHE"]["CE1"] = 2.1000;
	radii["PHE"]["CE2"] = 2.1000;
	radii["PHE"]["CZ"] = 2.1000;
	radii["PHE"]["C"] = 2.1000;
	radii["PHE"]["O"] = 1.6000;
	radii["PRO"]["N"] = 1.6000;
	radii["PRO"]["CD"] = 2.2350;
	radii["PRO"]["CA"] = 2.3650;
	radii["PRO"]["CB"] = 2.2350;
	radii["PRO"]["CG"] = 2.2350;
	radii["PRO"]["C"] = 2.1000;
	radii["PRO"]["O"] = 1.6000;
	radii["SER"]["N"] = 1.6000;
	radii["SER"]["H"] = 0.8000;
	radii["SER"]["CA"] = 2.3650;
	radii["SER"]["CB"] = 2.2350;
	radii["SER"]["OG"] = 1.6000;
	radii["SER"]["HG"] = 0.8000;
	radii["SER"]["C"] = 2.1000;
	radii["SER"]["O"] = 1.6000;
	radii["THR"]["N"] = 1.6000;
	radii["THR"]["H"] = 0.8000;
	radii["THR"]["CA"] = 2.3650;
	radii["THR"]["CB"] = 2.3650;
	radii["THR"]["OG1"] = 1.6000;
	radii["THR"]["HG1"] = 0.8000;
	radii["THR"]["CG2"] = 2.1650;
	radii["THR"]["C"] = 2.1000;
	radii["THR"]["O"] = 1.6000;
	radii["TRP"]["N"] = 1.6000;
	radii["TRP"]["H"] = 0.8000;
	radii["TRP"]["CA"] = 2.3650;
	radii["TRP"]["CB"] = 2.2350;
	radii["TRP"]["CG"] = 2.1000;
	radii["TRP"]["CD2"] = 2.1000;
	radii["TRP"]["CE2"] = 2.1000;
	radii["TRP"]["CE3"] = 2.1000;
	radii["TRP"]["CD1"] = 2.1000;
	radii["TRP"]["NE1"] = 1.6000;
	radii["TRP"]["HE1"] = 0.8000;
	radii["TRP"]["CZ2"] = 2.1000;
	radii["TRP"]["CZ3"] = 2.1000;
	radii["TRP"]["CH2"] = 2.1000;
	radii["TRP"]["C"] = 2.1000;
	radii["TRP"]["O"] = 1.6000;
	radii["TYR"]["N"] = 1.6000;
	radii["TYR"]["H"] = 0.8000;
	radii["TYR"]["CA"] = 2.3650;
	radii["TYR"]["CB"] = 2.2350;
	radii["TYR"]["CG"] = 2.1000;
	radii["TYR"]["CD1"] = 2.1000;
	radii["TYR"]["CE1"] = 2.1000;
	radii["TYR"]["CD2"] = 2.1000;
	radii["TYR"]["CE2"] = 2.1000;
	radii["TYR"]["CZ"] = 2.1000;
	radii["TYR"]["OH"] = 1.6000;
	radii["TYR"]["HH"] = 0.8000;
	radii["TYR"]["C"] = 2.1000;
	radii["TYR"]["O"] = 1.6000;
	radii["VAL"]["N"] = 1.6000;
	radii["VAL"]["H"] = 0.8000;
	radii["VAL"]["CA"] = 2.3650;
	radii["VAL"]["CB"] = 2.3650;
	radii["VAL"]["CG1"] = 2.1650;
	radii["VAL"]["CG2"] = 2.1650;
	radii["VAL"]["C"] = 2.1000;
	radii["VAL"]["O"] = 1.6000;
	radii["ACE"]["CH3"] = 2.1650;
	radii["ACE"]["C"] = 2.1000;
	radii["ACE"]["O"] = 1.6000;
	radii["AMN"]["CL"] = 2.1650;
	radii["AMN"]["C"] = 2.1000;
	radii["AMN"]["O"] = 1.6000;
	radii["CBX"]["N"] = 1.6000;
	radii["CBX"]["H"] = 0.8000;
	radii["CBX"]["CA"] = 2.1650;
	radii["O2"]["O1"] = 0.0000;
	radii["O2"]["O2"] = 0.0000;
	radii["CO"]["C"] = 0.0000;
	radii["CO"]["O"] = 0.0000;
	radii["ETH"]["C2"] = 2.1650;
	radii["ETH"]["C1"] = 2.2350;
	radii["ETH"]["O"] = 1.6000;
	radii["ETH"]["H"] = 0.8000;
	radii["ST2"]["ST2"] = 0.0000;
	radii["ST2"]["OX2"] = 0.0000;
	radii["ST2"]["HX1"] = 0.8000;
	radii["ST2"]["HX2"] = 0.8000;
	radii["ST2"]["LX1"] = 0.0000;
	radii["ST2"]["LX2"] = 0.0000;
	radii["OH2"]["OH2"] = 0.0000;
	radii["OH2"]["H1"] = 0.8000;
	radii["OH2"]["H2"] = 0.8000;
	radii["SO4"]["S"] = 1.8900;
	radii["SO4"]["O1"] = 1.6000;
	radii["SO4"]["O2"] = 1.6000;
	radii["SO4"]["O3"] = 1.6000;
	radii["SO4"]["O4"] = 1.6000;
	radii["COH"]["C"] = 2.1650;
	radii["COH"]["O"] = 1.6000;
	radii["COH"]["H"] = 0.8000;
	radii["T3P"]["OH2"] = 1.6000;
	radii["T3P"]["H1"] = 0.8000;
	radii["T3P"]["H2"] = 0.8000;
  return true;
}

// different
double resVDWRadii::maxSumRadii() {
	double m = 0.0;
	for (auto it1 = radii.begin(); it1 != radii.end(); it1++) {
		for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
			if(it2->second > m) {
				m = it2->second;
			}
		}
	}
	return m * 2;
}

// different
double resVDWRadii::getRadii(const string& resName, const string& atomName) {
	MstUtils::assert(resVDWRadii::radii.count(resName) > 0, "Residue type " + resName + " not found in resVDWRadii.");
	string aname = atomName;
	while (aname.length() >= 1) {
		if (resVDWRadii::radii[resName].count(aname) > 0) {
			if (atomName != aname) {
				cout << "WARNING: Replacing " + atomName + " with " + aname + ".\n";
			}
			return resVDWRadii::radii[resName][atomName];
		} else {
			aname = aname.substr(0, aname.size() - 1);
		}

	}
	MstUtils::assert(false, "Atom type " + atomName + " not found in resVDWRadii.");
  return double(0); //will never reach this
}

// different
double resVDWRadii::getRadii(Atom& a) {
	return getRadii(a.getParent()->getName(), a.getName());
}

double resVDWRadii::sumRadii(Atom& a1, Atom& a2) {
	return resVDWRadii::getRadii(a1) + resVDWRadii::getRadii(a2);
}

bool resVDWRadii::clash(Atom& a1, Atom& a2, double lb) {
	return (a1.distance(a2) < resVDWRadii::sumRadii(a1, a2) * lb);
}

bool resVDWRadii::contact(Atom& a1, Atom& a2, double lb, double ub) {
	double dist = a1.distance(a2);
	double s = resVDWRadii::sumRadii(a1, a2);
	return (dist >= s * lb && dist < s * ub);
}

bool resVDWRadii::independent(Atom& a1, Atom& a2, double ub) {
	return (a1.distance(a2) >= resVDWRadii::sumRadii(a1, a2) * ub);
}

atomInteraction resVDWRadii::interactionType(Atom& a1, Atom& a2, double lb, double ub) {
	double dist = a1.distance(a2);
	double s = resVDWRadii::sumRadii(a1, a2);
	if (dist < s * lb) {
		return CLASH;
	} else if (dist < s * ub) {
		return CONTACT;
	} else {
		return INDEPENDENT;
	}
}
