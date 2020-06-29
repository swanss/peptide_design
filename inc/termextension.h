//
//  TermExtension.h
//  TPD_swans
//
//  Created by Sebastian Swanson on 2/1/19.
//

#ifndef termextension_h
#define termextension_h

//mst dependencies
#include "msttypes.h"
#include "mstsystem.h"
#include "mstcondeg.h" //for ConFind
#include "dtermen.h" //for getContactsWith()
#include "mstfasst.h"
#include "mstmagic.h"

//tpd dependencies
#include "utilities.h" //for generateAllCombinationsKRes
#include "coverage.h" //to directly assess coverage
#include "secondarystructure.h" //to classify the secondary structure of the res
#include "structure_iter.h"
#include "vdwRadii.h"


/* TO DO:
 */
using namespace MST;

//// forward declarations
class TermExtension;
class seedTERM;
struct seed;

struct Seed {
    Structure* extended_fragment;
    mstreal rmsd; //rmsd for this specific match
    bool sequence_match; //true if match has same residue as the query at the central residue
    string secondary_structure_classification;
};


class seedTERM {
    friend class TermExtension;
public:
    /* Constructor */
    seedTERM();
    seedTERM(TermExtension* Fragmenter, vector<Residue*> all_res, bool search, bool seq_const);
    /*
     To do: fix copying Seeds. Currently, the pointers are copied, so if the parent object is deleted
     the pointers in the child object would be null.
     */
    seedTERM(const seedTERM& F);
    ~seedTERM();
    
    /* The seed type specifies how the residues from match proteins are selected to create the seed
     
     Many contact - each contact to the central residue that has flanking residues that overlap
     with another contact is combined to make a longer segment. Contacts that are not within this
     distance are saved as distinct segments, but the same seed.
     
     Peptide extension - only contacts that have flanking residues that overlap to the segment of the
     fragment in the original target are saved as seeds. */
    enum seedType {MANY_CONTACT, PEPTIDE_EXTENSION};
    
    // Fragment methods
    vector<Structure*> extendMatch(seedType seed_type, const Structure* match_structure, vector<int> match_idx, vector<vector<int>> seed_segments,vector<string> seed_sec_struct, mstreal RMSD, bool same_res);
    void writeExtendedFragmentstoPDB(string outDir, fstream& info_out, fstream& secstruct_out);
    void writeExtendedFragmentstoBIN(fstream& info_out, fstream& secstruct_out, StructuresBinaryFile* bin);
    //  void writeSecStruct(string outDir, fstream& secstruct_out);
    void deleteExtendedFragments();
    
    // Getters (objects)
    vector<int> getResIdx() const {return allResIdx;}
    const Structure& getStructure() const {return fragmentStructure;}
    vector<Residue*> getInteractingRes() const {return interactingRes;}
    fasstSolutionSet& getMatches() {return matches;}
    fasstSolution& getMatch(int i) {return matches[i];}
    //  vector<Structure*> getExtendedFragments() const {return extended_fragments;}
    
    
    // Getters (properties)
    string getCenResName() {return cenResName;}
    int getNumRes() {return fragmentStructure.residueSize();}
    int getNumMatches() {return matches.size();}
    int getSegmentNum(mstreal maxPeptideBondLen = 2.0);
    vector<int> getSegmentLens(mstreal maxPeptideBondLen = 2.0); //ordered by the residue indices
    vector<int> getResIdx() {return allResIdx;}
    Residue* getCenRes() {return cenRes;}
    int getCenResIdx() const {return cenResIdx;}
    TermExtension* getParent() const {return parent;}
    string getName() const {return name;};
    bool checkExtended() {return (!seeds.empty());}
    bool searchMode() {return search;}
    int getExtendedFragmentNum() {return num_extended_fragments;}
    mstreal getComplexity() {return effective_degrees_of_freedom;}
    mstreal getAdjRMSD() {return RMSD_cutoff;}
    mstreal getConstructionTime() {return construction_time;}
    
protected:
    // Called during Fragment construction
    void setName();
    
    // Called during Fragment extension
    void addResToStructure(vector<int> res_idx, bool keep_chain, const Structure* source_structure, Structure& recipient_structure);
    
    /* The extended fragment name is defined with the following components:
     -fragment name
     -seed flanking number
     -RMSD of the fragment to the match
     -seed number (unique ID)
     e.g. fragname_2-0.52319-132
     */
    string extendedFragmentName(mstreal RMSD);
    
private:
    TermExtension* parent;
    Structure fragmentStructure;
    Residue* cenRes;
    string cenResName;
    vector<Residue*> interactingRes; // These point to residues in the target structure
    vector<int> allResIdx;
    int cenResIdx; //index of the index of the central residue
    fasstSolutionSet matches;
    bool search;
    
    //  vector<Structure*> extended_fragments;
    //  vector<string> secondary_structure_classification;
    vector<Seed> seeds;
    mstreal RMSD_cutoff;
    mstreal effective_degrees_of_freedom;
    mstreal rmsd_adjust_factor;
    int num_extended_fragments;
    
    
    /* The fragment name is defined with the following components:
     -parent structure name
     -central residue chain AND number
     -chain AND number of all residues contacting the central residue
     -protein flanking residue number
     
     Ex: 5UUL-A12-A14A18A76-2
     */
    string name;
    mstreal construction_time;
    
    
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
    enum fragType {CEN_RES, CONTACT, ALL_COMBINATIONS, MATCH_NUM_REQ, COMPLEXITY_SCAN};
    
    /* Constructor */
    TermExtension(string fasstDBPath, string rotLib, vector<Residue*> selection, bool _verbose = true);
    
    /* Destructor */
    ~TermExtension();
    
    /* Fragmenter setters */
    void setTargetResidues(vector<Residue*> Selection);
    void setMatchReq(int matchReq) {match_req = matchReq;}
    void setFlankingNumber(int _flank) {flanking_res = _flank;}
    void setMinSeedLength(int len) {minimum_seed_length = len;}
    void setMaxRMSD(mstreal _max_rmsd) {max_rmsd = _max_rmsd;}
    void setSeqConst(bool _seq_const) {seq_const = _seq_const;}
    
    void resetFragments() {
        for (seedTERM* f : all_fragments) delete f;
        all_fragments.clear();
    }
    
    /* Fragmenter getters */
    Structure* getTarget() {return target_structure;}
    bool getAdaptiveRMSD() {return adaptive_rmsd;}
    int getFlank() const {return flanking_res;}
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
    
    void generateFragments(fragType option, bool search = true);
    
    /* Fragment getters */
    vector<seedTERM*> getFragments() {return all_fragments;}
    vector<Structure> getFragmentStructures();
    
    /* Fragment storage */
    // Note: if extendFragments has not been called, the seed number will not be stored
    void writeFragmentPDBs(string outDir);
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
    
    int extendFragment(seedTERM* f, seedTERM::seedType option, StructuresBinaryFile* bin, fstream& info, fstream& sec_struct);
    
    // Lower memory requirement
    void extendFragmentsandWriteStructures(seedTERM::seedType seed_type, string outDir = "");
    
    /* ExtendedFragment getters */
    vector<Structure*> getExtendedFragments();
    
    
protected:
    // Called by the constructor
    void set_params(string fasstDBPath, string rotLib);
    
    
    /* Called by extendFragments()
     Returns the indices of the residues that contact the central residue in the fasst target structure
     Each vector is a set of residues that will form a single seed segment
     NOTE: this is fundamentally different than the match_idx, because it is only the *contacting* residues */
    vector<int> identifySeedResidueIdx(const seedTERM* frag, const Structure* match_structure, vector<int> match_idx, int fasst_target_index);
    
    vector<int> getNonClashingResidueIdx(vector<int> seed_res_idx, const Structure* match_structure);
    
    vector<vector<int>> getSeedSegments(vector<int> seed_res_idx, const Structure* match_structure);
    
    vector<string> classifySegmentsInMatchProtein(Structure* S, vector<vector<int>> seed_segments);
    
private:
    // Main variables and parameters
    Structure* target_structure; // FULL structure (may include peptide chains)
    string structure_name;
    bool verbose;
    
    RotamerLibrary RL;
    string RL_path;
    
    FASST F;
    fasstSearchOptions foptsBase; // base FASST options that will be used with every search
    string fasstdbPath;
    int match_req; //if set to <= 0, does not influence process
    
    vector<Residue*> bindingSiteRes;
    mstreal cd_threshold, int_threshold, bbInteraction_cutoff;
    
    // Fragmentation parameters
    mstreal max_rmsd;
    int flanking_res; //must remain constant between fragment creation and seed extraction
    bool adaptive_rmsd, seq_const;
    
    // Fragments
    vector<seedTERM*> all_fragments;
    
    // Fragment extension parameters
    // variables stored for filtering seeds with clashes during TERM extension
    Structure target_BB_structure;
    AtomPointerVector target_BB_atoms;
    ProximitySearch* target_PS;
    vector<Structure*> target_structures;
    
    int seed_flanking_res, minimum_seed_length; // seed variables
    int extendedFragmentNumber; //counter
    
    secondaryStructureClassifier sec_structure;
    
};


#endif /* termextension_h */
