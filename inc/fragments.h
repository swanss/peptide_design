// 6/15/2020: This file is transferred and adapted from the structgen repo.

#ifndef FRAGMENTS_H
#define FRAGMENTS_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <climits>
#include <cfloat>

#include "msttypes.h"
#include "msttransforms.h"
#include "mstrotlib.h"
#include "mstcondeg.h"
#include "mstfasst.h"
#include "mstlinalg.h"

#include "utilities.h"
#include "params.h"

using namespace std;
using namespace MST;


class poseScore {

public:

	poseScore(int _search, int _adjust, bool _pass, double _rmsdAlign, double _rmsdOrig, double _rmsdWorst, double _rmsdAverage) {
		search = _search;
		adjust = _adjust;
		pass = _pass;
		rmsdAlign = _rmsdAlign;
		rmsdOrig = _rmsdOrig;
		rmsdWorst = _rmsdWorst;
		rmsdAverage = _rmsdAverage;
	}

	int getSearch() { return search; }
	int getAdjust() { return adjust; }
	int getPass() { return pass; }
	double getRMSDAlign() { return rmsdAlign; }
	double getRMSDOrig() { return rmsdOrig; }
	double getRMSDWorst() { return rmsdWorst; }
	double getRMSDAverage() { return rmsdAverage; }

	void setSearch(int _search) { search = _search; }
	void setAdjust(int _adjust) { adjust = _adjust; }
	void setPass(int _pass) { pass = _pass; }
	void setRMSDAdjust(double _rmsdAlign) { rmsdAlign = _rmsdAlign; }
	void setRMSDOrig(double _rmsdOrig) { rmsdOrig = _rmsdOrig; }
	void setRMSDWorst(double _rmsdWorst) { rmsdWorst = _rmsdWorst; }
	void setRMSDAverage(double _rmsdAverage) { rmsdAverage = _rmsdAverage; }

	bool operator < (const poseScore& ps) {
    if (search < ps.search) {
			return true;
		} else if (search > ps.search) {
			return false;
		}

		if (adjust < ps.adjust) {
			return true;
		} else if (adjust > ps.adjust) {
			return false;
		}

		if (pass > ps.pass) { // 1 > 0
			return true;
		} else if (pass < ps.pass) { // 0 < 1
			return false;
		}

		if (rmsdAlign < ps.rmsdAlign) {
			return true;
		} else if (rmsdAlign > ps.rmsdAlign) {
			return false;
		}

		if (rmsdAverage < ps.rmsdAverage) {
				return true;
		} else if (rmsdAverage > ps.rmsdAverage) {
				return false;
		}

		if (rmsdOrig < ps.rmsdOrig) {
			return true;
		} else {
			return false;
		}

	}

	void write(ostream& ss) {
		ss << "poseScore\n";
		ss << " search: " << search << endl;
		ss << " adjust: " << adjust << endl;
		ss << " pass: " << pass << endl;
		ss << " rmsdAlign: " << rmsdAlign << endl;
		ss << " rmsdOrig: " << rmsdOrig << endl;
		ss << " rmsdWorst: " << rmsdWorst << endl;
		ss << " rmsdAverage: " << rmsdAverage << endl;
	}

private:
	int search;
	int adjust;
	bool pass;
	double rmsdAlign; // rmsd of adjustment
	double rmsdOrig; // rmsd to original
	double rmsdWorst;
	double rmsdAverage;

};


class Fragment {
public:
	Fragment() { clear(); }

	// vector of residues to expand may be better because order of segments would be preserved
	// could have an order option
	// duplicates of residues in vector would not be kept
	Fragment(vector<Residue*>& _id, set<Residue*>& _contacts, vector<Residue*>& parentResidues, int numFlank, bool partial) {
		explode(_id, _contacts, parentResidues, numFlank, partial);
	}

	Fragment(vector<Residue*>& _id, set<Residue*>& _contacts, set<Residue*>& secondaryContacts, vector<Residue*>& parentResidues, int minNumFlank, int maxNumFlank, bool partial) {
		explodeVariable(_id, _contacts, secondaryContacts, parentResidues, minNumFlank, maxNumFlank, partial);
	}

	void clear () {
		id.clear();
		contacts.clear();
		residues.clear();
		expansions.clear();
	}

	void explode(vector<Residue*>& _id, set<Residue*>& _contacts, vector<Residue*>& parentResidues, int numFlank, bool partial);

	void explodeVariable(vector<Residue*>& _id, set<Residue*>& _contacts, set<Residue*>& secondaryContacts, vector<Residue*>& parentResidues, int minNumFlank, int maxNumFlank, bool partial);

	vector<vector<Residue*>> getSegments();

	bool isValid() {
		return (contacts.size() > 0 && residues.size() > 0);
	}

	bool operator < (const Fragment& other) const {
		return (contacts < other.contacts); // option for id, or residues
	}

	vector<Residue*> getID() {
		return id;
	}

	set<Residue*> getContacts() {
		return contacts;
	}

	map<Residue*, vector<Residue*>> getExpansions() {
		return expansions;
	}

	vector<Residue*> getExpansion(Residue* res) {
		MstUtils::assert(expansions.count(res) > 0, "Residue is not a key in expansions.");
		return expansions[res];
	}

	vector<Residue*> getResidues() {
		return residues;
	}

	// AtomPointerVector getAtoms() {
	// 	return apv;
	// }

	map<Residue*, int> residueIndexMapping() {
		map<Residue*, int> ret;
		for (int i = 0; i < residues.size(); i++) {
			ret[residues[i]] = i;
		}
		return ret;
	}

	vector<int> parentResidueIndices() {
		vector<int> ret;
		for (int i = 0; i < residues.size(); i++) {
			ret[i] = residues[i]->getResidueIndex();
		}
		return ret;
	}

	int fragmentResidueIndex(Residue* res) {
		for (int i = 0; i < residues.size(); i++) {
			if (res == residues[i]) {
				return i;
			}
		}
		MstUtils::assert(true, "Residue not found in fragment.");
        return -1;
	}

	vector<int> getNTerminiPositions();
	vector<int> getCTerminiPositions();
	int numSegments();

	Structure getStructure() {
		Structure s(residues);
		s.setName(toString());
		return s;
	}

	int numResidues() {
		return residues.size();
	}

	int numContacts() {
		return contacts.size();
	}

	friend ostream& operator <<(ostream& os, const Fragment& frag) {
		for (int i = 0; i < frag.id.size(); i++) {
			os << residueID(*frag.id[i]);
			if (i < frag.id.size() - 1) {
				os << "_";
			}
		}
        return os;
	}

	string toString() {
		stringstream ss;
		ss << *this;
		return ss.str();
	}

	void writePDB(const string& outFile) {
		Structure s = getStructure();
		s.writePDB(outFile);
	}

	void writeInfo(ostream& os) {
		os << toString() << endl;
		vector<Residue*> sorted = sortByResidueIndex(contacts);
		for (int i = 0; i < sorted.size(); i++) {
			os << "	" << residueID(*sorted[i]) << endl;
		}
	}

	void writeInfo(const string& outFile) {
		fstream of;
		MstUtils::openFile(of, outFile, fstream::out, "infoOut");
		writeInfo(of);
		of.close();
	}

private:
	vector<Residue*> id;
	set<Residue*> contacts;
	vector<Residue*> residues;
	map<Residue*, vector<Residue*>> expansions;
};


class Fragmenter {
	friend class Fragment;
public:
	Fragmenter() { reset(); }
	//Fragmenter(Fragment& _fragParams);
	Fragmenter(FragmentParams& fragParams);
	Fragmenter(const vector<Residue*>& _residues);
	Fragmenter(const vector<Residue*>& _residues, FragmentParams& _fragParams);
	Fragmenter (Structure& s);
	Fragmenter (Structure& s, FragmentParams& _fragParams);
	void init(const vector<Residue*>& _residues);
	void init(const vector<Residue*>& _residues, FragmentParams& _fragParams);
	void setResidues(const vector<Residue*>& _residues);
	void setResidues(Structure& s);
	void reset();
	void clear();
	void setParams(FragmentParams& _fragParams) { clear(); fragParams = _fragParams; }
	FragmentParams getParams() { return fragParams; }

	vector<Structure> getFragmentStructures();
	set<Fragment> getFragments() { return fragments; }
	map<Residue*, set<Residue*>> getContactMap() { return contactMap; }
	void buildContactMap(contactList& cl);
	void fragment(contactList& cl);

	friend ostream& operator <<(ostream& os, const Fragmenter& frags) {
		for (auto it = frags.fragments.begin(); it != frags.fragments.end(); it++) {
			os << *it << endl;
		}
		return os;
	}

	void writeContactMap(ostream& os, bool showNonContact = true);

private:
	vector<Residue*> residues;
	map<Residue*, set<Residue*>> contactMap;
	set<Fragment> fragments;
	FragmentParams fragParams;
};

struct fragmentAndSearchOutParams {
	bool protPDBOut, protContOut, dirListOut, infoOut, pdbOut, matchOut, seqOut, structOut, covOut, covResisOut, covBin, verbose;

	fragmentAndSearchOutParams() {
		protPDBOut = protContOut = dirListOut = infoOut = pdbOut = matchOut = seqOut = structOut = covOut = covResisOut = covBin = verbose = false;
	}

	fragmentAndSearchOutParams(int argc, char *argv[], MstOptions& opts, string ext = "") {
		opts.addOption("protPDBOut" + ext, "Write the parent protein(s) to the output directory(ies).", false);
		opts.addOption("protContOut" + ext, "Write the parent contacts(s) to the output directory(ies)", false);
		opts.addOption("pdbOut" + ext, "Write out the fragment PDB.", false);
		opts.addOption("infoOut"  + ext, "Write out the fragment information file.", false);
		opts.addOption("matchOut" + ext, "Write out the match file for each fragment.", false);
		opts.addOption("structOut" + ext, "Write out the match structures for each fragment.", false);
		opts.addOption("seqOut" + ext, "Write out the sequence file for each fragment.", false);
		opts.addOption("fragmentList" + ext, "Output list of directories for each PDB being fragmented.", false);
		opts.addOption("verbose", "Verbose option.", false);
		opts.addOption("covOut", "Write out elements covered by each fragment.");
		opts.addOption("covBin", "Write out elements covered to a binary file.");
		opts.addOption("covResisOut", "Write out the file mapping the residue indices of the matches to the covered elements.");
		opts.setOptions(argc, argv);

		protPDBOut = opts.isGiven("protPDBOut" + ext);
		protContOut = opts.isGiven("protContOut" + ext);
		pdbOut = opts.isGiven("pdbOut" + ext);
		infoOut = opts.isGiven("infoOut" + ext);
		matchOut = opts.isGiven("matchOut" + ext);
		structOut = opts.isGiven("structOut" + ext);
		seqOut = opts.isGiven("seqOut" + ext);
		verbose = opts.isGiven("verbose" + ext);
		dirListOut = opts.isGiven("fragmentList" + ext);
		covOut = opts.isGiven("covOut" + ext);
		covBin = opts.isGiven("covBin" + ext);
		covResisOut = opts.isGiven("covResisOut" + ext);
	}

};

contactList buildContacts(Structure& s, vector<Residue*>& residuesToInclude, RotamerLibrary* rl, contactParams& params, bool intra = false);

contactList buildContacts(Structure& s, RotamerLibrary* rl, contactParams& params, bool intra = false);

void splitContacts(Structure& s, vector<Residue*>& residuesToInclude, RotamerLibrary* rl, contactParams& params, bool intra, contactList& bbConts, contactList& bsConts, contactList&  sbConts, contactList& ssConts);

contactList contactListSubtract(contactList& c1, contactList& c2, bool directional = false);

contactList contactListSubtract(contactList& c, set<Residue*>& s, bool first, bool second);

#endif
