#ifndef termextension_h
#define termextension_h

//mst dependencies
#include "dtermen.h" //for getContactsWith()
#include "mstcondeg.h" //for ConFind
#include "mstfasst.h"
#include "mstmagic.h"
#include "mstsequence.h"
#include "mstsystem.h"
#include "msttypes.h"

//tpd dependencies
#include "benchmarkutilities.h"
#include "coverage.h" //to directly assess coverage
#include "freesasaext.h"
#include "secondarystructure.h" //to classify the secondary structure of the res
#include "structure_iter.h"
#include "vdwRadii.h"
#include "utilities.h" //for generateAllCombinationsKRes

/* TO DO:
 */
using namespace MST;

// forward declarations
struct seed;
class seedTERM;
class TermExtension;

class TEParams;
class interfaceCoverage;

class nearbySeedScore;

struct Seed {
    Structure* extended_fragment;
    int match_number; //the rank of this match
    mstreal rmsd; //rmsd for this specific match
    bool sequence_match; //true if match has same residue as the query at the central residue
    string secondary_structure_classification;
    int num_fragment_matches;
};


class seedTERM {
    friend class TermExtension;
public:
    /* Constructor */
    seedTERM();
    seedTERM(TermExtension* Fragmenter, vector<Residue*>_interactingRes, vector<int> _flankResPerSegment, bool search, mstreal max_rmsd, bool discard_match_data = false);
//    /*
//     To do: fix copying Seeds. Currently, the pointers are copied, so if the parent object is deleted
//     the pointers in the child object would be null.
//     */
//    seedTERM(const seedTERM& F);
    ~seedTERM();
    
    /* The seed type specifies how the residues from match proteins are selected to create the seed
     
     Many contact - each contact to the central residue that has flanking residues that overlap
     with another contact is combined to make a longer segment. Contacts that are not within this
     distance are saved as distinct segments, but the same seed.ß
     
     Peptide extension - only contacts that have flanking residues that overlap to the segment of the
     fragment in the original target are saved as seeds. */
    enum seedType {MANY_CONTACT, PEPTIDE_EXTENSION};
    
    // Fragment methods
    vector<Structure*> extendMatch(seedType seed_type, const Structure* match_structure, vector<int> match_idx, vector<vector<int>> seed_segments,vector<string> seed_sec_struct, int match_number, mstreal RMSD, int target_ID, bool same_res);
    void writeExtendedFragmentstoPDB(string outDir, fstream& info_out, fstream& secstruct_out, interfaceCoverage* IC);
    void writeExtendedFragmentstoBIN(fstream& info_out, fstream& secstruct_out, StructuresBinaryFile* bin, interfaceCoverage* IC);
    //  void writeSecStruct(string outDir, fstream& secstruct_out);
    void deleteExtendedFragments();
    
    bool operator<(seedTERM* other) {
        /**
        Two step comparison:
        1. Compare number of residues in fragment structure
        2. If 1 is tied, compare number of matches.
         */
        if (fragmentStructure.residueSize() == other->fragmentStructure.residueSize()) {
            if ((!search) || (!other->search)) MstUtils::error("Cannot compare two seedTERMs with equal number of residues unless they were also searched for matched","seedTERM::operator<");
            return tD.numMatches() < other->tD.numMatches();
        } else {
            return fragmentStructure.residueSize() < other->fragmentStructure.residueSize();
        }
    }
    
    bool operator==(seedTERM* other) {
        /**
         Two fragments are the same only if they have all the same residues
         */
        if (fragmentStructure.residueSize() != other->fragmentStructure.residueSize()) return false;
        else {
            for (int i = 0; i < fragmentStructure.residueSize(); i++) {
                if (&fragmentStructure.getResidue(i) != &other->fragmentStructure.getResidue(i)) return false;
            }
        }
        return true;
    }
    
    // Getters (objects)
    vector<int> getResIdx() const {return allResIdx;}
    const Structure& getStructure() const {return fragmentStructure;}
    vector<Residue*> getInteractingRes() const {return interactingRes;}
    vector<int> getFlankingResPerInteractingRes() const {return flankResPerInteractingRes;}
    
    // Potentially unsafe to call if "discard_match_data" is set to true
    fasstSolutionSet& getMatches() {return tD.getMatches();}
    fasstSolution& getMatch(int i) {return tD.getMatch(i);}
    
    
    // Getters (properties)
    string getCenResName() {return cenResName;}
    int getNumRes() {return fragmentStructure.residueSize();}
    int getNumMatches() {return numMatches;}
    int getNumMatchesPreFiltering() {return numMatchesPreHomologyFiltering;}
    int getSegmentNum(mstreal maxPeptideBondLen = 2.0);
    vector<int> getSegmentLens(mstreal maxPeptideBondLen = 2.0); //ordered by the residue indices
    vector<int> getResIdx() {return allResIdx;}
    Residue* getCenRes() const {return cenRes;}
    int getCenResIdx() const {return cenResIdx;}
    TermExtension* getParent() const {return parent;}
    string getName() const {return name;};
    bool checkExtended() {return (!seeds.empty());}
    bool searchMode() {return search;}
    int getExtendedFragmentNum() {return num_extended_fragments;}
    int getSeedResNum() {return num_seed_residues;}
    mstreal getComplexity() {return effective_degrees_of_freedom;}
    mstreal getMaxRMSD() {return max_RMSD_cutoff;}
    mstreal getAdaptiveRMSD() {return adaptive_RMSD_cutoff;}
    mstreal getSearchRMSD() {return RMSD_cutoff;}
    mstreal getConstructionTime() {return construction_time;}
    
protected:
    // Called during Fragment construction
    void setName();
    
    // Called during Fragment extension
    void addResToStructure(vector<int> res_idx, bool keep_chain, const Structure* source_structure, Structure& recipient_structure);
    
    /* The extended fragment name is defined with the following components:
     -fragment name
     -match protein ID in db
     -seed flanking number
     -RMSD of the fragment to the match
     -seed number (unique ID)
     e.g. fragname_2-0.52319-132
     */
    string extendedFragmentName(mstreal RMSD, int match_protein_ID);
    
private:
    TermExtension* parent;
    Structure fragmentStructure;
    Residue* cenRes;
    string cenResName;
    vector<Residue*> interactingRes; // These point to residues in the target structure
    vector<int> flankResPerInteractingRes;
    vector<int> allResIdx;
    int cenResIdx; //index of the central residue in fragmentStructure
//    string seed_chain_ID = "0";
    
    int numMatchesPreHomologyFiltering;
    int numMatches;
    termData tD;
    
    bool search;
    
    //  vector<Structure*> extended_fragments;
    //  vector<string> secondary_structure_classification;
    vector<Seed> seeds;
    mstreal max_RMSD_cutoff;
    mstreal adaptive_RMSD_cutoff;
    mstreal RMSD_cutoff; //this is the RMSD cutoff that is actually applied in the search
    mstreal effective_degrees_of_freedom;
    mstreal rmsd_adjust_factor;
    int num_extended_fragments;
    int num_seed_residues;
    
    
    /* The fragment name is defined with the following components:
     -parent structure name
     -central residue chain AND number
     -chain AND number of all residues contacting the central residue
     -the number of flanking residues for central res/contacting residues (no delim)
     
     Ex: 5UUL-A12-A14A18A76-3223
     */
    string name;
    mstreal construction_time;
};

class TEParams {
    friend seedTERM;
    friend TermExtension;
public:
    TEParams(string params_file_path) {
        setParamsFromFile(params_file_path);
    };
    
    void setParamsFromFile(string params_file_path);
    
    void printValues();
    
    string getConfigFile() {return config_file;}

    enum fragType {CEN_RES, ALL_COMBINATIONS, ADAPTIVE_SIZE, ADAPTIVE_LENGTH, ADAPTIVE_LENGTH_FIXED_RMSD, COMPLEXITY_SCAN};
    
    /**
     Sequence constraints are used to pre-filter which positions will be searched in the database. Only those segments that match the
     constraints will be considered in the structural search
     
     NONE - no sequence constraint is applied
     CEN_RES_ONLY - the central residue must have the same amino acid as the query in the fragment
     ALL_RES - all residues of the match must have the same amino acid as the query
     BLOSUM62 - all residues with scores of 0 or greater (meaning they are neutral or better) are considered equivalent. The match must
     have all residues equivalent to those at the query positions
     HYBRID - same as BLOSUM62, but the central residue must be an exact match.
     */
    
    enum seqConstType {NONE, CEN_RES_ONLY, ALL_RES, BLOSUM62, HYBRID};
    
    fragType fragment_type = ADAPTIVE_LENGTH;
    mstreal cd_threshold = .01;
    mstreal cdSeqConst_threshold = .01;
    mstreal int_threshold = .01;
    mstreal bbInteraction_cutoff = 3.5;
    mstreal max_rmsd = 1.2;
    int flanking_res = 2;
    int match_req = -1; //only considered if value is > 0
    bool adaptive_rmsd = false;
    seqConstType seq_const = NONE;
    string config_file = "";
    int seed_flanking_res = 2;
    mstreal homology_cutoff = 0.5;
    bool allow_sidechain_clash = true;
    mstreal freedom_cutoff = 0.4;
    mstreal relSASA_cutoff = -1.0;
    bool verbose = false;
    string acceptable_aa_substitutions_path = "";
    
    map<string,vector<string>> acceptable_aa_substitutions;
    
protected:
    void loadAASubMap();
};


class TermExtension {
    friend class seedTERM;
    
    /* The fragmenter stores the objects/parameters used to define fragments from a structure of interest.
     The fragments are searched against the provided database to identify matches and these matches are
     used to identify new residues that could be used for design, or "seeds".
     */
public:
    /*  FRAGMENT_TYPE DEFINITIONS
     
     A fragment is specified by a set of interacting residues. Each of these residues, plus the union
     of their +/- flanking residues, are included in the final fragment structure. Note that regardless
     of the specific type, every fragment will always include the binding site residue upon which
     it is centered.
     
     Below, each type specifies how the interacting residues are chosen. In each type, the union
     of the flanking residues is *always* included.
     
     Types:
     Central residue - just the central residue.
     
     Contact - all residues that contact the binding site residue (i.e. its neighbourhood).
     
     In the following fragment types, some combination of the neighbourhood residues are chosen.
     
     All combinations - Every possible combination of neighbourhood residues, R, is used to make
     a fragment. Results in a fragment with all the residues in R, R choose R-1 fragments with R-1 interacting
     residues, ... , R choose 0 residues with no interacting residues (just the binding site residue).
     SUPER MEMORY/STORAGE INTENSIVE
     
     Match Number Requirement - The set of interacting residues that gives the largest possible
     fragment (defined by *total* number of residues) while still having at least the required
     number of matches, N. (Or at least a fragment close to that fragment). Results in one fragment
     per binding site residue.
     
     Complexity Scan - same underlying logic as Match Number Requirement, with the following
     two differences: 1) a fragment is constructed each time an interacting residue is added and
     2) residues are added until no more remain. This will result in R number fragments, where R is
     the number of residues in the neighbourhood.
     
     */
//    enum fragType {CEN_RES, ALL_COMBINATIONS, ADAPTIVE_SIZE, ADAPTIVE_LENGTH, COMPLEXITY_SCAN};
    
    /* Constructor */
    TermExtension(string fasstDBPath, string rotLib, vector<Residue*> selection, TEParams& _params);
    
    /* Destructor */
    ~TermExtension();
    
    /* Fragmenter setters */
    void setTargetResidues(vector<Residue*> Selection);
    void setMatchReq(int matchReq) {params.match_req = matchReq;}
    void setFlankingNumber(int _flank) {params.flanking_res = _flank;}
    void setMinSeedLength(int len) {minimum_seed_length = len;}
    void setMaxRMSD(mstreal _max_rmsd) {params.max_rmsd = _max_rmsd;}
    void setAdaptiveRMSD(bool _adaptive_rmsd) {params.adaptive_rmsd = _adaptive_rmsd;}
//    void setSeqConst(bool _seq_const) {params.seq_const = _seq_const;}
    void setIC(interfaceCoverage* _IC) {IC = _IC;}
    void setWriteAllFiles(string _extFragDir, string _wholeMatchProteinDir, string _seedDir) {
        storeAll = true;
        extFragDir = _extFragDir;
        wholeMatchProteinDir = _wholeMatchProteinDir;
        seedDir = _seedDir;
    }
    void setVCalcRadius(mstreal qR) {vCalc.setQueryRadius(qR);}
  
    void resetFragments() {
        for (seedTERM* f : all_fragments) delete f;
        all_fragments.clear();
    }
    
    /* Fragmenter getters */
    Structure* getTarget() {return target_structure;}
    bool getAdaptiveRMSD() {return params.adaptive_rmsd;}
    int getFlank() const {return params.flanking_res;}
    string getName() {return structure_name;}
    int getExtendedFragmentNumber() {return extendedFragmentNumber;}
    
    /* Store */
    void storeParameters(string Dir);
    void storePeptideProteinContacts(string Dir);
    
    /*  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -*/
    
    /* Methods for making fragments from the binding site residues
     
     option - The type of fragments that will be created.
     
     search - By default, the fragments are searched during their construction, but if that is not
     desired, it can be controlled by changing this variable to false. */
    
    void generateFragments(bool search = true);
    
    /* Fragment getters */
    vector<seedTERM*> getFragments() {return all_fragments;}
    vector<Structure> getFragmentStructures();
    
    /* Fragment storage */
    // Note: if extendFragments has not been called, the seed number will not be stored
    void writeFragmentPDBs(string outDir, string binnedDataPath = "");
    void writeFragmentClassification(string outDir);
    
    /*  -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -*/
    
    // Method for extending fragments to generate seeds
    /* In order to find new structures, the matches to each fragment will be examined. The residue
     of the match that corresponds to the central residue in the original fragment will be used
     to search for new contacts (that were not present in the original fragment). These contacts
     will be expanded on both ends and combined with other overlapping contacts or stored, each
     contact forming its own seed.
     
     Note: if store is selected, the extended fragments are written to a file and then deleted. This
     helps reduce the memory requirement.
     
     Options to build in:
     -Only new contacts but not all new residues (extension of a chain) (extension)
     Generating structures from the selected residues
     */
    //  void extendFragments(seedType seed_type, bool store = 0, string outDir = "", string file_type = "bin");
    
    int extendFragment(seedTERM* f, seedTERM::seedType option, StructuresBinaryFile* bin, fstream& matchInfo, fstream& info, fstream& sec_struct);
    
    // Lower memory requirement
    void extendFragmentsandWriteStructures(seedTERM::seedType seed_type, string outDir = "");
    
    /* ExtendedFragment getters */
    vector<Structure*> getExtendedFragments();
    
    string seed_chain_ID = "0";
    
protected:
    // Called by the constructor
    void set_params(string fasstDBPath, string rotLib);
    
    /**
     Takes a single seedTERM and finds all possible expansions. Expansions can be generated by 1) incrementing the number of flanking
     residues on a segment or 2) adding a residue that contacts the center one (which creates a new segment)
     
     @param fragToExpand A fragment that we want to generate seedTERM expansions of
     @param contactingRes Residues that could potentially be added to the fragment
     @param search Controls whether the fragments are searched against the DB when constructed. In some cases not all fragments need to be compared when selecting the best, so avoiding the search can save time
     @return All viable expansions of the fragToExpand
     */
    vector<seedTERM*> generateAllExpansions(seedTERM* fragToExpand, vector<Residue*> contactingRes, bool search);
    
    /**
     Takes seedTERMs (from generateAllExpansions) and selects the one with 1) the most residues and 2) as a tie-breaker, the most matches.
     If no seedTERM has sufficient matches, none are returned.
     
     @param fragExpansions The seedTERMs that will be compared
     @return the seedTERM that satisifies the criteria, or if there is none, nullptr
     */
    seedTERM* bestFragment(vector<seedTERM*> fragExpansions);
    
    /* Called by extendFragments()
     Returns the indices of the residues that contact the central residue in the fasst target structure
     Each vector is a set of residues that will form a single seed segment
     NOTE: this is fundamentally different than the match_idx, because it is only the *contacting* residues */
    vector<int> identifySeedResidueIdx(const seedTERM* frag, const Structure* match_structure, vector<int> match_idx, int fasst_target_index, string& cenResAA, fstream& match);
    
    vector<int> getNonClashingResidueIdx(vector<int> seed_res_idx, const Structure* match_structure, fstream& match);
    
    vector<vector<int>> getSeedSegments(vector<int> seed_res_idx, const Structure* match_structure);
    
    vector<string> classifySegmentsInMatchProtein(Structure* S, vector<vector<int>> seed_segments);
    
private:
    // Main variables and parameters
    Structure* target_structure; // FULL structure
    string structure_name;
    bool verbose;
    
    RotamerLibrary RL;
    string RL_path;
    
    FASST F;
    fasstSearchOptions foptsBase; // base FASST options that will be used with every search
    string fasstdbPath;
    
    interfaceCoverage* IC = nullptr; //if not null, will check if seed overlaps peptide when deciding whether to store
    
    vector<Residue*> bindingSiteRes;
    
    // Fragments
    vector<seedTERM*> all_fragments;
    
    // Fragment extension parameters
    // variables stored for filtering seeds with clashes during TERM extension
    // NOT necessarily all BB atoms: depends on user params
    Structure target_BB_structure;
    AtomPointerVector target_BB_atoms;
    ProximitySearch* target_PS;
    vector<vector<Atom*>> target_structures_atoms;
    
    int minimum_seed_length; // seed variables
    int extendedFragmentNumber; //counter
    
    void setAAToSeqContProp();
    map<string,string> aaToSeqContProp; //for getting the keys to the fasst property
    
    secondaryStructureClassifier sec_structure;
    
    TEParams params; //handles all parameters
    
    bool storeAll = false; //controls whether extra files are written out
    string extFragDir = "";
    string wholeMatchProteinDir = "";
    string seedDir = "";
    
    sasaCalculator sasaCalc;
    map<Residue*,mstreal> res2relSASA;
  
    volumeCalculator vCalc;
};

class nearbySeedScore {
public:
    nearbySeedScore(vector<Residue*> _selectedRes, string seedBinaryPath, string _binnedDataPath, int _maxNumMatches, mstreal _cR, mstreal vR) : selectedRes(_selectedRes), seeds(seedBinaryPath), binnedDataPath(_binnedDataPath), maxNumMatches(_maxNumMatches), cR(_cR) {
        if (selectedRes.empty()) MstUtils::error("Must provide at least one residue","nearbySeedScore::nearbySeedScore");
        
        target = selectedRes[0]->getStructure();
        for (Residue* R : selectedRes) for (Atom* A : R->getAtoms()) if (A->getName() == "CA") selCAatoms.push_back(A);
        if (selCAatoms.size() != selectedRes.size()) MstUtils::error("Mismatch number of CA atoms and residues in the provided selection");
        psCA = ProximitySearch(selCAatoms,cR/2);
        
        vCalc = volumeCalculator(target,vR);
        
        if (binnedDataPath != "") {
            binnedData.readHistFile(binnedDataPath);
            binnedData.exceptOutOfRangeQueries();
        }
        
        for (Residue* R : selectedRes) selRes2NearbySeedResCount[R] = 0.0;
        }
    
    void scoreResidues();
    
    void writeResidueInfo(string prefix, set<Residue*> bindingSiteRes);
protected:
private:
    mstreal maxNumMatches = 100;
    mstreal cR = 12.0; //count radius
    
    vector<Residue*> selectedRes; //residues to be scored
    Structure* target; //whole target structure
    vector<Atom*> selCAatoms;
    ProximitySearch psCA;
    
    StructureIterator seeds;
    
    volumeCalculator vCalc;
    
    string binnedDataPath;
    oneDimBinnedData binnedData;
    
    map<Residue*,mstreal> selRes2NearbySeedResCount;
    map<Residue*,mstreal> selRes2FracUnoccupiedVol;
    map<Residue*,mstreal> selRes2Score;
};


#endif /* termextension_h */
