// 6/15/2020: This file is transferred from the structgen repo.

#ifndef PARAMS_H
#define PARAMS_H

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <climits>

#include "mstoptions.h"
#include "mstfasst.h"
#include "msttypes.h"

#include "utilities.h"

using namespace std;
using namespace MST;


RotamerLibrary* getRotamerLibrary(int argc, char *argv[]);


struct rmsdParams {
public:
	double maxRMSD;
	double persLen;
	int minGapConn;

	rmsdParams(double _maxRMSD = 1.0, double _persLen = 15, int _minGapConn = 1) {
		maxRMSD = _maxRMSD;
		persLen = _persLen;
		minGapConn = _minGapConn;
	}

	rmsdParams(int argc, char *argv[], MstOptions& opts, string ext = "", double _maxRMSD = 1.0, double _persLen = 15, int _minGap = 1) {
		opts.addOption("maxRMSD" + ext, "Maximum RMSD for " + ext, false);
		opts.addOption("persLen" + ext, "Persistence length for " + ext, false);
		opts.addOption("minGap" + ext, "Minium gap length to consider as separate segment for " + ext, false);
		opts.setOptions(argc, argv);
		maxRMSD = opts.getReal("maxRMSD" + ext, _maxRMSD);
		persLen = opts.getReal("persLen" + ext, _persLen);
		minGapConn = opts.getInt("minGap" + ext, _minGap);
	}

	double rmsdCutoff(vector<int>& segLens) {
		return RMSDCalculator::rmsdCutoff(segLens, maxRMSD, persLen);
	}

	double rmsdCutoff(Structure& s) {
		// option not to split??
		return RMSDCalculator::rmsdCutoff(s.reassignChainsByConnectivity(), maxRMSD, persLen);
	}

	double rmsdCutoff(vector<Residue*>& queryResidues, int minLen = 3) {
		Structure s(queryResidues);
		s = s.reassignChainsByConnectivity();
		vector<int> lens = chainLengths(s);
		for (int i = 0; i < lens.size(); i++) {
			lens[i] = min(lens[i], minLen);
		}
		return RMSDCalculator::rmsdCutoff(lens, maxRMSD, persLen);
	}

	double rmsdCutoff(vector<Residue*>& queryResidues, Structure& s) {
		vector<int> segLens = effectiveSegLengths(s, queryResidues, minGapConn);
		return RMSDCalculator::rmsdCutoff(segLens, maxRMSD, persLen);
	}

	void write(ostream& ss) {
// 		ss << "rmsdParams\n";
		ss << "maxRMSD " << maxRMSD << endl;
		ss << "persLen " << persLen << endl;
		ss << "minGapConn " << minGapConn << endl;
	}

};


struct contactParams {
	double bbDeg;
	double sbDeg;
	double ssDeg;

	contactParams(double _bbDeg = 3.5, double _sbDeg = 0.03, double _ssDeg = 0.03) {
		bbDeg = _bbDeg;
		sbDeg = _sbDeg;
		ssDeg = _ssDeg;
	}

	contactParams(int argc, char *argv[], MstOptions& opts, string ext = "", double _bbDeg = 3.5, double _sbDeg = 0.02, double _ssDeg = 0.02) {
		opts.addOption("bbDeg" + ext, "", false);
		opts.addOption("sbDeg" + ext, "", false);
		opts.addOption("ssDeg" + ext, "", false);
		opts.setOptions(argc, argv);
		bbDeg = opts.getReal("bbDeg" + ext, _bbDeg);
		sbDeg = opts.getReal("sbDeg" + ext, _sbDeg);
		ssDeg = opts.getReal("ssDeg" + ext, _ssDeg);
	}

	void write(ostream& ss) {
// 		ss << "contactParams\n";
		ss << "bbDeg " << bbDeg << endl;
		ss << "sbDeg " << sbDeg << endl;
		ss << "ssDeg " << ssDeg << endl;
	}
};


struct FragmentParams {
	int minNumFlank;
	int maxNumFlank;
	bool partialFlank;
	bool pair; // should be removed  - separate pair and residue fragmenters
	bool mustContact;
	bool variableFlank;
	// pairs ?
	// self ?

	FragmentParams() { init(2, false, false, false, false); }

	FragmentParams(int _numFlank, bool _pair, bool _partialFlank = false, bool _mustContact = false, bool _variableFlank = false) {
		init(_numFlank, _numFlank, _pair, _partialFlank, _mustContact, _variableFlank);
	}

	FragmentParams(int _minNumFlank, int _maxNumFlank, bool _pair, bool _partialFlank = false, bool _mustContact = false, bool _variableFlank = false) {
		init(_minNumFlank, _maxNumFlank, _partialFlank, _mustContact, _variableFlank);
	}

	void init(int _minNumFlank, int _maxNumFlank, bool _pair, bool _partialFlank = false, bool _mustContact = false, bool _variableFlank = false) {
		minNumFlank = _minNumFlank;
		maxNumFlank = _maxNumFlank;
		partialFlank = _partialFlank;
		pair = _pair;
		mustContact = _mustContact;
		variableFlank = _variableFlank;
	}

	FragmentParams(int argc, char *argv[], MstOptions& opts, string ext = "") {
		opts.addOption("minNumFlank" + ext, "", false);
		opts.addOption("maxNumFlank" + ext, "", false);
		opts.addOption("numFlank" + ext, "", false);
		opts.addOption("pair" + ext, "", false);
		opts.addOption("partialFlank" + ext, "", false);
		opts.addOption("mustContact" + ext, "", false);
		opts.addOption("variableFlank" + ext, "", false);

		MstUtils::assert(((opts.isGiven("minNumFlank" + ext) || opts.isGiven("maxNumFlank" + ext)) != opts.isGiven("numFlank" + ext)), "minNumFlank or maxNumFlank cannot be given with numFlank.");
		if (opts.isGiven("numFlank" + ext)) {
			minNumFlank = opts.getInt("numFlank" + ext, 1);
			maxNumFlank = opts.getInt("numFlank"  + ext, 1);
		} else {
			minNumFlank = opts.getInt("minNumFlank" + ext, 1);
			maxNumFlank = opts.getInt("maxNumFlank" + ext, 1);
		}
		pair = opts.isGiven("pair");
		mustContact = opts.isGiven("mustContact");
		partialFlank = opts.isGiven("partialFlank" + ext); // allow shorter fragments at segment ends
		variableFlank = opts.isGiven("variableFlank" + ext);
	}

	void write(ostream& ss) {
		ss << "minNumFlank " << minNumFlank << endl;
		ss << "maxNumFlank " << maxNumFlank << endl;
		ss << "partialFlank " << partialFlank << endl;
	}
};


struct DesParams2 {
	rmsdParams rParamsSS;
	rmsdParams rParamsSB;
	rmsdParams rParamsBB;
	contactParams contParams;
	int numFlank;
	int maxNumMatches;
	bool terminal;

	DesParams2() {}

	DesParams2(int argc, char *argv[], MstOptions& opts, string ext = "") {
		opts.addOption("numFlank" + ext, "", false);
		opts.addOption("terminal" + ext, "", false);
		opts.addOption("maxNumMatches" + ext, "", false);
		opts.setOptions(argc, argv);
		numFlank = opts.getInt("numFlank" + ext, 1);
		maxNumMatches = opts.getInt("maxNumMatches" + ext, 5000);
		terminal = !opts.isGiven("noEnds" + ext);
		contParams = contactParams(argc, argv, opts, ext);
		rParamsBB = rmsdParams(argc, argv, opts, ext + "BB", 0.80, 15.0, 1);
		rParamsSB = rmsdParams(argc, argv, opts, ext + "SB", 0.90, 15.0, 1);
		rParamsSS = rmsdParams(argc, argv, opts, ext + "SS", 0.90, 15.0, 1);
	}

	void write(ostream& ss) {
// 		ss << "desParams\n";
		ss << "numFlank " << numFlank << endl;
		ss << "maxNumMatches " << maxNumMatches << endl;
		ss << "terminal " << terminal << endl;
		contParams.write(ss);
		ss << "BB RMSD params:\n";
		rParamsBB.write(ss);
		ss << "SB RMSD params:\n";
		rParamsSB.write(ss);
		ss << "SS RMSD params:\n";
		rParamsSS.write(ss);
	}
};


// change to search params
struct fasstParams {
	double maxRMSD;
	double persLen;
	double maxSeqIdent;
	int sufficientNumMatches;
	int minNumMatches;
	int maxNumMatches;
// 	int minGap;
// 	int maxGap;
	bool bbRMSD;

	int minGapConn;

	fasstParams() {};

	fasstParams(double _maxRMSD, double _persLen, int _minGapConn = 1, double _maxSeqIdent = 1.0, int _sufficientNumMatches = INT_MAX, int _minNumMatches = 0, int _maxNumMatches = INT_MAX, bool _bbRMSD = true) {
		maxRMSD = _maxRMSD;
		persLen = _persLen;
		minGapConn = _minGapConn;
		maxSeqIdent = _maxSeqIdent;
		sufficientNumMatches = _sufficientNumMatches;
		minNumMatches = _minNumMatches;
		maxNumMatches = _maxNumMatches;
// 		minGap = _minGap;
// 		maxGap = _maxGap;
		bbRMSD = _bbRMSD;
	}

	fasstParams(int argc, char *argv[], MstOptions& opts, string ext = "") {
		opts.addOption("maxRMSD" + ext, "", true);
		opts.addOption("persLen" + ext, "", true);
		opts.addOption("maxSeqIdent" + ext, "", false);
		opts.addOption("sufficientNumMatches" + ext, "", false);
		opts.addOption("minNumMatches" + ext, "", false);
		opts.addOption("maxNumMatches" + ext, "", false);
// 		opts.addOption("minGap" + ext, "", false);
// 		opts.addOption("maxGap" + ext, "", false);
		opts.addOption("caRMSD" + ext, "", false);

		opts.setOptions(argc, argv);
		maxRMSD = opts.getReal("maxRMSD" + ext);
		persLen = opts.getReal("persLen" + ext);
		maxSeqIdent = opts.getReal("maxSeqIdent" + ext, 100) / 100.0;
		sufficientNumMatches = opts.getInt("sufficientNumMatches" + ext, INT_MAX);
		minNumMatches = opts.getInt("minNumMatches" + ext, 0);
		maxNumMatches = opts.getInt("maxNumMatches" + ext, INT_MAX);
// 		minGap = opts.getInt("minGap" + ext, 0);
// 		maxGap = opts.getInt("maxGap" + ext, INT_MAX);
		bbRMSD = !opts.isGiven("caRMSD" + ext);
	}

	void setFasst(FASST& fasst) {
		fasst.setRedundancyCut(maxSeqIdent);
		fasst.setSufficientNumMatches(sufficientNumMatches);
		fasst.setMinNumMatches(minNumMatches);
		fasst.setMaxNumMatches(maxNumMatches);
// 		fasst.setMinGap(minGap);
// 		fasst.setMaxGap(maxGap);
		if (bbRMSD) {
			fasst.setSearchType(FASST::FULLBB);
		} else {
			fasst.setSearchType(FASST::CA);
		}
	}

	void setFasst(FASST& fasst, Structure& query) {
		setFasst(fasst);
		setRMSD(fasst, query);
	}

	void setFasst(FASST& fasst, Structure& parent, vector<Residue*> queryResidues) {
		setFasst(fasst);
		setRMSD(fasst, parent, queryResidues);
	}

	void setFasst(FASST& fasst, vector<int>& segLens) {
		setFasst(fasst);
		setRMSD(fasst, segLens);
	}

	void setRMSD(FASST& fasst, Structure& query) {
		Structure tmp = query.reassignChainsByConnectivity();
		double rmsdCutoff = RMSDCalculator::rmsdCutoff(tmp, maxRMSD, persLen);
		fasst.setRMSDCutoff(rmsdCutoff);
	}

	void setRMSD(FASST& fasst, Structure& s, vector<Residue*> queryResidues) {
		vector<int> segLens = effectiveSegLengths(s, queryResidues, minGapConn);
		setRMSD(fasst, segLens);
	}

	void setRMSD(FASST& fasst, vector<int>& segLens) {
		double rmsdCutoff = RMSDCalculator::rmsdCutoff(segLens, maxRMSD, persLen);
		fasst.setRMSDCutoff(rmsdCutoff);
	}

	void write(ostream& ss) {
		ss << "maxRMSD " << maxRMSD << endl;
		ss << "persLen " << persLen << endl;
		ss << "minGapConn " << minGapConn << endl;
		ss << "maxSeqIdent " << maxSeqIdent << endl;
		ss << "sufficientNumMatches " << sufficientNumMatches << endl;
		ss << "minNumMatches " << minNumMatches << endl;
		ss << "maxNumMatches " << maxNumMatches << endl;
// 		ss << "minGap " << minGap << endl;
// 		ss << "maxGap " << maxGap << endl;
		ss << "bbRMSD " << bbRMSD << endl;
	}
};


// change to search params
struct MasterParams {
	double maxRMSD;
	double persLen;
	int topN;
	bool bbRMSD;
	double distCut; // get rid of

	MasterParams() {};

	MasterParams(double _maxRMSD, double _persLen, int _topN = INT_MAX, bool _bbRMSD = true, double _distCut = 40.0) {
		maxRMSD = _maxRMSD;
		persLen = _persLen;
		topN = _topN;
		bbRMSD = _bbRMSD;
		double distCut = _distCut;
	}

	MasterParams(int argc, char *argv[], MstOptions& opts, string ext = "") {
		opts.addOption("maxRMSD" + ext, "", true);
		opts.addOption("persLen" + ext, "", true);
		opts.addOption("topN" + ext, "", false);
		opts.addOption("caRMSD" + ext, "", false);
		opts.addOption("distCut" + ext, "", false);

		opts.setOptions(argc, argv);
		maxRMSD = opts.getReal("maxRMSD" + ext);
		persLen = opts.getReal("persLen" + ext);
		topN = opts.getInt("topN" + ext, INT_MAX);
		bbRMSD = !opts.isGiven("caRMSD" + ext);
		distCut = opts.getReal("distCut" + ext, 40.0);
	}

	void write(ostream& ss) {
		ss << "maxRMSD " << maxRMSD << endl;
		ss << "persLen " << persLen << endl;
		ss << "topN " << topN << endl;
		ss << "bbRMSD " << bbRMSD << endl;
		ss << "distCut " << distCut << endl;
	}
};


struct DesParams {
	double interCut;
	int numFlank;
	bool terminal;
	double maxRMSD;
	double persLen;
	bool bbRMSD;
	int maxNumMatches;
	double maxPepBond;
	double distCut;
	int minSep;

	DesParams() { //double _interCut, int _numFlank = 1, double _maxRMSD = 0.8, double _persLen = 15, int _maxNumMatches = 1, bool _terminal = true, bool _bbRMSD = true, double _maxPepBond = 2.0, double _distCut = 40.0, int minSep = 3) {
		interCut = -1;
		numFlank = 1;
		terminal = true;
		maxRMSD = 0.7;
		persLen = 15;
		bbRMSD = true;
		maxNumMatches = 1;
		maxPepBond = 2.0;
	}

	bool isSet() {
		return (interCut > 0);
	}


	DesParams(int argc, char *argv[], MstOptions& opts, string ext = "") {
		opts.addOption("interCut" + ext, "", true);
		opts.addOption("numFlank" + ext, "", false);
		opts.addOption("noEnds" + ext, "", false);
		opts.addOption("maxRMSD" + ext, "", false);
		opts.addOption("persLen" + ext, "", false);
		opts.addOption("caRMSD" + ext, "", false);
		opts.addOption("maxNumMatches" + ext, "", false);
		opts.addOption("maxPepBond" + ext, "", false);
		opts.addOption("minSep" + ext, "", false);

		opts.setOptions(argc, argv);
		interCut = opts.getReal("interCut" + ext);
		numFlank = opts.getInt("numFlank" + ext, 1);
		terminal = !opts.isGiven("noEnds" + ext);
		maxRMSD = opts.getReal("maxRMSD" + ext, 0.8);
		persLen = opts.getReal("persLen"  + ext, 10);
		bbRMSD = !opts.isGiven("caRMSD" + ext);
		maxNumMatches = opts.getInt("maxNumMatches" + ext, 1);
		maxPepBond = opts.getReal("maxPepBond" + ext, 2.0);
		minSep = opts.getInt("minSep", 3);
	}

	void write(ostream& ss) {
		ss << "interCut " << interCut << endl;
		ss << "terminal " << terminal << endl;
		ss << "numFlank " << numFlank << endl;
		ss << "maxRMSD " << maxRMSD << endl;
		ss << "persLen " << persLen << endl;
		ss << "bbRMSD " << bbRMSD << endl;
		ss << "maxNumMatches " << maxNumMatches << endl;
		ss << "maxPepBond " << maxPepBond << endl;
		ss << "minSep " << minSep << endl;
	}
};


// struct LoopParams {
// 	int minLoopLen;
// 	int maxLoopLen;
// 	int minNumLoops;
// 	int maxNumLoops;
//
// 	int matchLen;
// 	MasterParams masterParams;
// 	vector<int> topN; // the topN parameter given to master in the searches for termini matches
// 	vector<int> nextTopN; // the topN extensions to go on to the next round - will be sorted
// 	vector<int> nextGapLens;
// 	double ident; // the maximum identity allowed in the matches
// 	double maxRMSDRed;
// 	double persLenRed;
//
// 	int numFlankRed;
//
// 	double vdwUnfused;
// 	double vdwFused;
//
// 	int extLen;
// 	int minExtLen;
// 	int maxExtLen;
// 	int maxNumExt;
//
// 	double devDist;
// 	double devDistDec;
//
// 	bool intra;
// 	rmsdParams rParams;
// 	DesParams2 desParams;
// 	int scoreMethod;
// 	bool localDes;
// 	vector<int> structChains;
//
// 	int loopLenStep;
// 	string badLoopDir;
//
// 	LoopParams() {}
//
// 	LoopParams(int _minLoopLen, int _maxLoopLen, rmsdParams& _rParams, DesParams2& _desParams,
// 				int _minNumLoops = 1, int _maxNumLoops = INT_MAX, vector<int> _topN = {INT_MAX}, vector<int> _nextTopN = {50}, double _ident = 0.4, int _numFlankRed = 11,
// 				int _matchLen = 2, int _extLen = 3, int _minExtLen = 2, int _maxNumExt = 1, double _devDist = 4.0, double _devDistDec = 2.0,
// 				vector<int> _nextGapLens = vector<int>(), bool _intra = false, int _scoreMethod = 1, bool _localDes = false, double _vdwUnfused = 0.75, double _vdwFused = 0.80) {
//
// 		minLoopLen = _minLoopLen;
// 		maxLoopLen = _maxLoopLen;
// 		minNumLoops = _minNumLoops; // will keep searching until this many have been found - if over this many have been found - only keep this many
// 		maxNumLoops = _maxNumLoops; // this is the maximum number of loops to store // should be store in some sort of priority queue
//
// 		rParams = _rParams;
// 		desParams = _desParams;
//
// 		ident = _ident;
// 		numFlankRed = _numFlankRed;
// 		matchLen = _matchLen;
// 		extLen = _extLen;
// 		minExtLen = _minExtLen;
//
// 		maxNumExt = _maxNumExt;
// 		devDist = _devDist;
// 		devDistDec = _devDistDec;
// 		scoreMethod = _scoreMethod;
// 		localDes = _localDes;
// 		vdwUnfused = _vdwUnfused;
// 		vdwFused = _vdwFused;
//
// 		intra = _intra;
// 		assignVectorParams(_topN, _nextTopN, _nextGapLens);
// 	}
//
//
// 	void assignVectorParams(vector<int>& _topN, vector<int>& _nextTopN, vector<int>& _nextGapLens) {
// 		topN = _topN;
// 		if (topN.size() < maxNumExt) {
// 			int lastIdx = topN.size() - 1;
// 			int last = topN[lastIdx];
// 			for (int i = lastIdx + 1; i < maxNumExt; i++) {
// 				topN.push_back(last);
// 			}
// 		}
//
// 		nextTopN = _nextTopN;
// 		if (nextTopN.size() < maxNumExt - 1) {
// 			int lastIdx = nextTopN.size() - 1;
// 			int last = nextTopN[lastIdx];
// 			for (int i = lastIdx + 1; i < maxNumExt - 1; i++) {
// 				nextTopN.push_back(last);
// 			}
// 		}
//
// 		if (_nextGapLens.size() == 0) {
// 			for (int i = 1; i < maxNumExt; i++) {
// 				nextGapLens.push_back(INT_MAX);
// 			}
// 		} else {
// 			nextGapLens = _nextGapLens;
// 			if (nextGapLens.size() < maxNumExt - 1) {
// 				int lastIdx = nextGapLens.size() - 1;
// 				int last = nextGapLens[lastIdx];
// 				for (int i = lastIdx + 1; i < maxNumExt - 1; i++) {
// 					nextGapLens.push_back(last);
// 				}
// 			}
// 		}
// 	}
//
//
// 	LoopParams(int argc, char *argv[], MstOptions& opts, string ext = "") {
// 		opts.addOption("pdbs", "Database of PDB structures to be used in MASTER searches. The default is /home/ironfs/scratch/grigoryanlab/cmack2357/bc-30-clean/list_pdb.", false);
//
// 		opts.addOption("minLoopLen", "Minimum loop length. This is the minimum number of new residues that can be added to the structure. The default is 1", false);
// 		opts.addOption("maxLoopLen", "Maximum loop length. This is the maximum number of new residues that can be added to the structure. The default is 8", false);
// 		opts.addOption("matchLen", "The match length (number of resideus in the overlap region) of the extension. The default is 2.", false);
// 		opts.addOption("extLen", "The number of residues to be added beyond the overlap region. The default is 2.", false);
// 		opts.addOption("minExtLen", "The minimum number of residues to be added beyond the overlap region. Only extensions with at least this many residues will be considered. This will only be considered if for some reason the match connoted be extend the requested number of residues. The default is 1.", false);
// 		opts.addOption("maxExtLen", "The maximum number of residues to be added beyond the overlap region. Only extensions with at most this many residues will be considered. This will only be considered if for some reason the match connoted be extend the requested number of residues. The default is 3.", false);
//
// 		opts.addOption("topN", "Consider at most this many matches for each loop extension. Can be a single number used at all levels of recursion or a comma separated list of numbers (e.g. 100,500,100). The default is 1000", false);
// 		opts.addOption("nextTopN", "At most this many matches should go on to next round. Can be a single number used at all levels of recursion or a comma separated list of numbers (e.g. 100,500,100). The default is ", false);
// 		opts.addOption("ident", "Percent identity for redundancy removal when search for loops. The default is 40%.", false);
// 		opts.addOption("maxRMSDRed", "Maximum RMSD used in structural redundancy removal. The default is 0.5.", false);
// 		opts.addOption("persLenRed", "Persistence length used in structural redundancy removal. The default is 10.", false);
// 		opts.addOption("numFlankRed", "Number of flanking residues for sequence redundancy removal. The default is 11 residues on either side of the central residues.", false);
//
// 		opts.addOption("minNumLoops", "Minimum number of solutions (loops) to consider. If this many solutions have not been found, it will go into the next level of recursion (unless it has reached the maximum level of recursion).", false);
// 		opts.addOption("maxNumLoops", "Maximum number of solutions (loops) to consider. The program will quit when finds this many solutions.", false);
// 		opts.addOption("maxNumExt", "The maximum number of extensions to consider while building loops. The default is 1.", false);
//
// 		opts.addOption("devDist", "This determines how much further the C-N distance of the first (non-closing) extension can be compared to the original C-N distance. The default is 4.0A", false);
// 		opts.addOption("devDistDec", "This determines how much the distance cutoff between C- and N-termini should be decreased per extension. The default is 2.0.", false);
//
// 		opts.addOption("intra", "Base the final design score on all contacts with the loop region. Default will be just those contacts between the loop region and the existing region (not including contacts between residues in the loop region).", false);
// 		opts.addOption("scoreMethod", "Contact based structural designability: 0, Contact based sequence/structure based compatability: 1");
// 		opts.addOption("localDes", "Create backbone fragments and search for structural designability with this many flanking residues. If 0 (default) is given, no searches will be done.", false);
// 		opts.addOption("loopLenStep", "Maximum number of loops that are recorded at a given length. Will be visited in order of RMSD for each loop length.", false);
// 		opts.addOption("badLoopDir", "Directory to write out loops that pass all filters except for designability.", false);
// 		opts.addOption("vdwUnfused", "Scale factor for van der waals radii in clash checking for fused loop..", false);
// 		opts.addOption("vdwFused", "Scale factor for van der waals radii in clash checking for unfused loop.", false);
// 		opts.addOption("structChains", "Chains to leave out of sequence scoring.", false);
//
// 		opts.setOptions(argc, argv);
//
// 		string pdbs = opts.getString("pdbs", "/home/ironfs/scratch/grigoryanlab/cmack2357/bc-30-clean/list_pdb");
//
// 		minLoopLen = opts.getInt("minLoopLen", 1);
// 		maxLoopLen = opts.getInt("maxLoopLen", 8);
// 		matchLen = opts.getInt("matchLen", 2);
// 		extLen = opts.getInt("extLen", 2);
// 		minExtLen = opts.getInt("minExtLen", 2);
// 		maxExtLen = opts.getInt("maxExtLen", 3);
//
// 		ident = opts.getReal("ident", 50.0) / 100.0;
// 		maxRMSDRed = opts.getReal("maxRMSDRed", 0.5);
// 		persLenRed = opts.getReal("persLenRed", 10);
// 		numFlankRed = opts.getInt("numFlankRed", 11);
//
// 		devDist = opts.getReal("devDist", 4.0);
// 		devDistDec = opts.getReal("devDistDec", 2.0);
//
// 		minNumLoops = opts.getInt("minNumLoops", 1);
// 		maxNumLoops = opts.getInt("maxNumLoops", INT_MAX);
// 		maxNumExt = opts.getInt("maxNumExt", 1);
//
// 		intra = opts.isGiven("intra");
// 		scoreMethod = opts.getInt("scoreMethod", 1);
// 		localDes = opts.isGiven("localDes");
// 		loopLenStep = opts.getInt("loopLenStep", INT_MAX);
// 		badLoopDir = opts.getString("badLoopDir", "");
//
// 		vdwFused = opts.getReal("vdwFused", 0.8);
// 		vdwUnfused = opts.getReal("vdwUnfused" , 0.75);
//
// 		string structChainsString = opts.getString("structChains", "");
// 		vector<string> elems = MstUtils::split(structChainsString, ",");
// 		for (int i = 0; i < elems.size(); i++) {
// 			structChains.push_back(stoi(elems[i]) - 1);
// 		}
//
// 		string topNString = opts.getString("topN", "10000");
// 		elems = MstUtils::split(topNString, ",");
// 		vector<int> _topN;
// 		for (int i = 0; i < elems.size(); i++) {
// 			_topN.push_back(stoi(elems[i]));
// 		}
// 		// extend last value if less than the number of extensions
//
// 		string nextTopNString = opts.getString("nextTopN", "50");
// 		elems = MstUtils::split(nextTopNString, ",");
// 		vector<int> _nextTopN;
// 		for (int i = 0; i < elems.size(); i++) {
// 			_nextTopN.push_back(stoi(elems[i]));
// 		}
//
// 		string nextGapLensString = opts.getString("nextGapLens", "1000");
// 		elems = MstUtils::split(nextGapLensString, ",");
// 		vector<int> _nextGapLens;
// 		for (int i = 0; i < elems.size(); i++) {
// 			_nextGapLens.push_back(stoi(elems[i]));
// 		}
//
// 		assignVectorParams(_topN, _nextTopN, _nextGapLens);
//
// 		rParams = rmsdParams(argc, argv, opts);
// 		// use contactParams
// 		// use fragmentParams
// 		// use rmsdParams for anchor
// 		desParams = DesParams2(argc, argv, opts, "Des");
// 	}
//
// 	void write(ostream& ss) {
// 		ss << "LoopParams\n";
// 		ss << "minLoopLen " << minLoopLen << endl;
// 		ss << "maxLoopLen " << maxLoopLen << endl;
// 		ss << "matchLen " << matchLen << endl;
// 		ss << "extLen " << extLen << endl;
// 		ss << "minExtLen " << minExtLen << endl;
// 		ss << "devDist " << devDist << endl;
// 		ss << "devDistDec " << devDistDec << endl;
// 		ss << "minNumLoops " << minNumLoops << endl;
// 		ss << "maxNumLoops " << maxNumLoops << endl;
// 		ss << "maxNumExt " << maxNumExt << endl;
// 		ss << "intra " << intra << endl;
// 		ss << "ident " << ident << endl;
// 		ss << "maxRMSDRed " << maxRMSDRed << endl;
// 		ss << "persLenRed " << persLenRed << endl;
// 		ss << "structChains "; writeVector(ss, structChains); ss << endl;
// 		ss << "topN "; writeVector(ss, topN); ss << endl;
// 		ss << "nextTopN "; writeVector(ss, nextTopN); ss << endl;
// 		ss << "nextGapLens "; writeVector(ss, nextGapLens); ss << endl;
// 		ss << "scoreMethod " << scoreMethod << endl;
// 		ss << "localDes " << localDes << endl;
// 		ss << "loopLenStep " << loopLenStep << endl;
// 		ss << "vdwUnfused " << vdwUnfused << endl;
// 		ss << "vdwFused " << vdwFused << endl;
// 		ss << "********\n";
// 		ss << "\nDesParams\n";
// 		desParams.write(ss);
// 		ss << "********\n";
// 		ss << "********\n";
// 	}
// };

#endif
