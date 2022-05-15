#include "fragments.h"
#undef assert


void Fragment::explode(vector<Residue*>& _id, set<Residue*>& _contacts, vector<Residue*>& parentResidues, int numFlank, bool partialFlank) {
	id = _id;
	contacts = _contacts;
	residues.clear();
	expansions.clear();

	set<int> fragIndices;
	for (auto const& central: contacts) {
		int centralIndex = central->getResidueIndex();
		int start = max(0, centralIndex - numFlank);
		set<int> centralExpansion;
		for (int i = centralIndex - 1; i >= start; i--) {
			if (!isBonded(*parentResidues[i], *parentResidues[i + 1])) {
				break;
			}
			centralExpansion.insert(i);
		}

		if (centralExpansion.size() < numFlank && !partialFlank) { // under size AND not terminal
			fragIndices.clear();
			break;
		}

		centralExpansion.insert(centralIndex);
		int end = min(centralIndex + numFlank, (int) parentResidues.size() - 1);
		for (int i = centralIndex + 1; i <= end; i++) {
			if (!isBonded(*parentResidues[i - 1], *parentResidues[i])) {
				break;
			}
			centralExpansion.insert(i);
		}

		if (centralExpansion.size() < (2 * numFlank + 1) && !partialFlank) { // under size AND not terminal
			fragIndices.clear();
			break;
		}
		fragIndices.insert(centralExpansion.begin(), centralExpansion.end());
		for (auto const& it: centralExpansion) {
			expansions[central].push_back(parentResidues[it]);
		}
	}

	if (fragIndices.size() > 0) {
		for (auto const& it: fragIndices) {
			residues.push_back(parentResidues[it]);
		}
	}
}


void Fragment::explodeVariable(vector<Residue*>& _id, set<Residue*>& _contacts, set<Residue*>& secondaryContacts, vector<Residue*>& parentResidues, int minNumFlank, int maxNumFlank, bool partialFlank) {
	id = _id;
	contacts = _contacts;
	residues.clear();
	expansions.clear();
	set<int> fragIndices;
	for (auto const& central: contacts) {
		int centralIndex = central->getResidueIndex();
		set<int> centralExpansion;

		int maxStart = max(0, centralIndex - minNumFlank);
		int minStart = max(0, centralIndex - maxNumFlank);
		for (int i = centralIndex - 1; i >= minStart; i--) {
			if (!isBonded(*parentResidues[i], *parentResidues[i + 1]) ||
			 		(i < maxStart && secondaryContacts.count(parentResidues[i]) == 0)) {
				break;
			}
			centralExpansion.insert(i);
		}


		if (centralExpansion.size() < minNumFlank && !partialFlank) { // under size AND not partialFlank
			fragIndices.clear();
			break;
		}

		int minEnd = min((int) parentResidues.size() - 1, centralIndex + minNumFlank);
		int maxEnd = min((int) parentResidues.size() - 1, centralIndex + maxNumFlank);

		centralExpansion.insert(centralIndex);
		for (int i = centralIndex + 1; i <= maxEnd; i++) {
			if (!isBonded(*parentResidues[i - 1], *parentResidues[i]) ||
					(i > minEnd && secondaryContacts.count(parentResidues[i]) == 0)) {
				break;
			}
			centralExpansion.insert(i);
		}

		if (centralExpansion.size() < (2 * minNumFlank + 1) && !partialFlank) { // under size AND not terminal
			fragIndices.clear();
			break;
		}
		fragIndices.insert(centralExpansion.begin(), centralExpansion.end());
		for (auto const& res: centralExpansion) {
			expansions[central].push_back(parentResidues[res]);
		}
	}

	if (fragIndices.size() > 0) {
		for (auto const& it: fragIndices) {
			residues.push_back(parentResidues[it]);
		}
	}
}


vector<vector<Residue*>> Fragment::getSegments() {
	vector<vector<Residue*>> segResidues;
	vector<Residue*> tmp;
	for (int i = 0; i < residues.size(); i++) {
		tmp.push_back(residues[i]);
		if (i == residues.size() - 1 || !isBonded(*residues[i], *residues[i + 1])) {
			segResidues.push_back(tmp);
			tmp.clear();
		}
	}
	return segResidues;
}


int Fragment::numSegments() {
	int n = 0;
	for (int i = 0; i < residues.size(); i++) {
		if (i == 0 || !isBonded(*residues[i - 1], *residues[i])) {
			n++;
		}
	}
	return n;
}

vector<int> Fragment::getNTerminiPositions() {
	vector<int> result;
	for (int i = 0; i < residues.size(); i++) {
		if (i == 0 || !isBonded(*residues[i - 1], *residues[i])) {
			result.push_back(i);
		}
	}
	return result;
}

vector<int> Fragment::getCTerminiPositions() {
	vector<int> result;
	for (int i = 0; i < residues.size(); i++) {
		if (i == residues.size() - 1 || !isBonded(*residues[i], *residues[i + 1])) {
			result.push_back(i);
		}
	}
	return result;
}


void Fragmenter::init(const vector<Residue*>& _residues) {
	reset();
	setResidues(_residues);
}

void Fragmenter::init(const vector<Residue*>& _residues, FragmentParams& _fragParams) {
	init(_residues);
	fragParams = _fragParams;
}


Fragmenter::Fragmenter(FragmentParams& _fragParams) {
	reset();
	fragParams = _fragParams;
}

Fragmenter::Fragmenter(const vector<Residue*>& _residues) {
	init(_residues);
}

Fragmenter::Fragmenter(const vector<Residue*>& _residues, FragmentParams& _fragParams) {
	init(_residues, _fragParams);
}

Fragmenter::Fragmenter(Structure& s) {
	init(s.getResidues());
}

Fragmenter::Fragmenter(Structure& s, FragmentParams& _fragParams) {
	init(s.getResidues(), _fragParams);
}

void Fragmenter::setResidues(const vector<Residue*>& _residues) {
	clear();
	residues = _residues;
}

void Fragmenter::setResidues(Structure& s) {
	setResidues(s.getResidues());
}

void Fragmenter::reset() {
	fragParams.init(2, 2, false, false, true, false);
	clear();
}



void Fragmenter::clear() {
	contactMap.clear();
	fragments.clear();
}

vector<Structure> Fragmenter::getFragmentStructures() {
	vector<Structure> fragStructs;
	for (auto it = fragments.begin(); it != fragments.end(); it++) {
		Fragment frag = *it;
		fragStructs.push_back(frag.getStructure());
	}
	return fragStructs;
}

void Fragmenter::buildContactMap(contactList& cl) { // set ??
	clear();
	if (fragParams.pair) {
		for (int i = 0; i < cl.size(); i++) {
			Residue* resA = cl.residueA(i);
			Residue* resB = cl.residueB(i);
			if (resA->getResidueIndex() < resB->getResidueIndex()) {
				contactMap[resA].insert(resB);
			} else {
				contactMap[resB].insert(resA);
			}
		}
	} else {
		for (int i = 0; i < cl.size(); i++) {
			contactMap[cl.residueA(i)].insert(cl.residueB(i));
			contactMap[cl.residueB(i)].insert(cl.residueA(i));
		}
	}
}

void Fragmenter::fragment(contactList& cl) {
	buildContactMap(cl);
	set<Residue*> secondaryContacts;
	if (fragParams.pair) {
		for (auto it1 = contactMap.begin(); it1 != contactMap.end(); it1++) {
			Residue* res1 = it1->first;
			for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
				Residue* res2 = *it2;
				vector<Residue*> id = {res1, res2};
				id = sortByResidueIndex(id);
				set<Residue*> contacts = {res1, res2};
				if (fragParams.variableFlank) {
					secondaryContacts.insert(contactMap[res1].begin(), contactMap[res1].end());
					secondaryContacts.insert(contactMap[res2].begin(), contactMap[res2].end());
				}
				Fragment frag(id, contacts, secondaryContacts, residues, fragParams.minNumFlank, fragParams.maxNumFlank, fragParams.partialFlank);
				if (frag.isValid()) {
					fragments.insert(frag);
				}
			}
		}
	} else {
		for (int i = 0; i < residues.size(); i++) {
			Residue* res1 = residues[i];
			vector<Residue*> id = {res1};
			set<Residue*> contacts = {res1};
			if (contactMap.count(res1) > 0) {
				contacts.insert(contactMap[res1].begin(), contactMap[res1].end());
			}
			if (!fragParams.mustContact || contacts.size() > 1) {
				if (fragParams.variableFlank) {
					for (auto const& res: contacts) {
						secondaryContacts.insert(contactMap[res].begin(), contactMap[res].end());
					}
				}
				Fragment frag(id, contacts, secondaryContacts, residues, fragParams.minNumFlank, fragParams.maxNumFlank, fragParams.partialFlank);
				if (frag.isValid()) {
					fragments.insert(frag);
				}
			}
		}
	}
}

void Fragmenter::writeContactMap(ostream& os, bool showNonContact) {
	for (int i = 0; i < residues.size(); i++) {
		Residue* res = residues[i];
		if (showNonContact || contactMap.count(res) > 0) {
			os << residueID(*res) << " -> ";
		}
		if (contactMap.count(res) > 0) {
			set<Residue*>& contacts = contactMap[res];
			vector<Residue*> sorted = sortByResidueIndex(contacts);
			for (int j = 0; j < sorted.size(); j++) {
				os << residueID(*sorted[j]);
				if (j < sorted.size() - 1) {
					os << " ";
				}
			}
		}
		os << endl;
	}
}


contactList buildContacts(Structure& s, vector<Residue*>& residuesToInclude, RotamerLibrary* rl, contactParams& params, bool intra) {
	set<Residue*> toCheck(residuesToInclude.begin(), residuesToInclude.end());

	contactList allConts = backboneContacts(s, params.bbDeg);
	contactList contsBB;
	set<pair<Residue*, Residue*> > seenContacts;
	for (int i = 0; i < allConts.size(); i++) {
		Residue* resA = allConts.residueA(i);
		Residue* resB = allConts.residueB(i);
		int numIncluded = toCheck.count(resA) + toCheck.count(resB);
		if (numIncluded == 1 || (numIncluded == 2 && intra)) { // exaclty one residue must be from the loop region - could think about allowing two for intra
			contsBB.addContact(resA, resB, allConts.degree(i));
			seenContacts.insert(make_pair(resA, resB));
		}
	}

	ConFind confind(rl, s);
	contactList tmp1 = confind.getContacts(residuesToInclude, params.ssDeg);
	contactList tmp2 = confind.getInterference(residuesToInclude, params.sbDeg);
	allConts = contactListUnion({tmp1, tmp2});
	contactList contsS;
	for (int i = 0; i < allConts.size(); i++) {
		Residue* resA = allConts.residueA(i);
		Residue* resB = allConts.residueB(i);
		int numIncluded = toCheck.count(resA) + toCheck.count(resB);
		if (seenContacts.count(make_pair(resA, resB)) == 0 && (numIncluded == 1 || (numIncluded == 2 && intra))) {
			contsS.addContact(resA, resB, allConts.degree(i));
		}
	}
	return contactListUnion({contsBB, contsS});
}


contactList buildContacts(Structure& s, RotamerLibrary* rl, contactParams& params, bool intra) {
	vector<Residue*> residuesToInclude = s.getResidues();
	return buildContacts(s, residuesToInclude, rl, params, intra);
}


void splitContacts(Structure& s, vector<Residue*>& residuesToInclude, RotamerLibrary* rl, contactParams& params, bool intra, contactList& bbConts, contactList& bsConts, contactList& sbConts, contactList& ssConts, bool verbose) {
	set<Residue*> toCheck(residuesToInclude.begin(), residuesToInclude.end());

	bbConts = contactList(); // BB of toCheck and BB of s
	bsConts = contactList(); // BB of toCheck and SC of s
	sbConts = contactList(); // SC of toCheck and BB of s (unless intra, then SC of toCheck and BB of toCheck)
	ssConts = contactList(); // SC of toCheck and SC of s

	ConFind confind(rl, s);
	contactList conts = confind.getContacts(residuesToInclude, params.ssDeg);
    if (verbose) cout << "sidechain-sidechain contacts:" << conts.size() << endl;
	for (int i = 0; i < conts.size(); i++) {
		Residue* resA = conts.residueA(i);
		Residue* resB = conts.residueB(i);
        if (verbose) cout << *resA << " and " << *resB << endl;
		int numIncluded = toCheck.count(resA) + toCheck.count(resB);
		if ((!intra && numIncluded == 1) || (intra && numIncluded == 2)) {
            if (verbose) cout << "adding..." << endl;
			ssConts.addContact(resA, resB, conts.degree(i));
		}
	}

	conts = confind.getInterference(residuesToInclude, params.sbDeg);
    if (verbose) cout << "sidechain-backbone/backbone-sidechain contacts:" << conts.size() << endl;
	for (int i = 0; i < conts.size(); i++) {
		Residue* resA = conts.residueA(i); // side chain from A to backbone of B
		Residue* resB = conts.residueB(i);
        if (verbose) cout << *resA << " and " << *resB << endl;
        //the "source" residues is the one whose sidechain interferes with the others backbone
        bool sourceIsInToCheck = toCheck.count(resA) == 1; //the source is in the residues to be designed
        int numIncluded = toCheck.count(resA) + toCheck.count(resB);
        if (!sourceIsInToCheck && ((intra && numIncluded == 2) || (!intra && numIncluded == 1))) { // BB of toCheck
            if (verbose) cout << "adding... (and swapping direction)" << endl;
            bsConts.addContact(resB, resA, conts.degree(i)); // side chain coming from existing (target) region
        }
		else if (sourceIsInToCheck && ((intra && numIncluded == 2) || (!intra && numIncluded == 1))) {
			if (verbose) cout << "adding..." << endl;
            sbConts.addContact(resA, resB, conts.degree(i)); // side chain from toCheck
		}
	}
    
    conts = confind.getBBInteraction(residuesToInclude, params.bbDeg);
//    conts = backboneContacts(s, params.bbDeg);
    if (verbose) cout << "backbone-backbone contacts:" << conts.size() << endl;
	for (int i = 0; i < conts.size(); i++) {
		Residue* resA = conts.residueA(i);
		Residue* resB = conts.residueB(i);
        if (verbose) cout << *resA << " and " << *resB << endl;
		int numIncluded = toCheck.count(resA) + toCheck.count(resB);
		if ((!intra && numIncluded == 1) || (intra && numIncluded == 2)) { // exaclty one residue must be from the loop region - could think about allowing two for intra
			if (verbose) cout << "adding..." << endl;
            bbConts.addContact(resA, resB, conts.degree(i));
		}
	}
}

// remove elements of c2 from c1
contactList contactListSubtract(contactList& c1, contactList& c2, bool directional) {
	set<pair<Residue*, Residue*>> pairs2;
	for (int i = 0; i < c2.size(); i++) {
		pairs2.insert(make_pair(c2.residueA(i), c2.residueB(i)));
		if (!directional) {
			pairs2.insert(make_pair(c2.residueB(i), c2.residueA(i)));
		}
	}

	contactList diff;
	for (int i = 0; i < c1.size(); i++) {
		Residue* resA = c1.residueA(i);
		Residue* resB = c1.residueB(i);
		if (pairs2.count(make_pair(resA, resB)) == 0) {
				diff.addContact(resA, resB, c1.degree(i));
		}
	}
	return diff;
}


contactList contactListSubtract(contactList& c, set<Residue*>& s, bool first, bool second) {
	contactList diff;
	for (int i = 0; i < c.size(); i++) {
		if ((!first || s.count(c.residueA(i)) == 0) &&
				(!second || s.count(c.residueB(i)) == 0)) {
					diff.addContact(c.residueA(i), c.residueB(i), c.degree(i));
		}
	}
	return diff;
}
