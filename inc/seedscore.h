#ifndef seedscore_h
#define seedscore_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <set>
#include <unordered_map>
#include <fstream>
#include <math.h>
#define _USE_MATH_DEFINES

#include "msttypes.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "msttransforms.h"
#include "mstrotlib.h"
#include "mstfasst.h"
#include "dtermen.h"
#include "dtermenscorer.h"

#include "seedgraph.h"

#include "utilities.h"
#include "params.h"
#include "fragments.h"

using namespace std;
using namespace MST;


/**
 Creates a semicolon-separated string representing the scores for each residue
 in the score map.
 */
string writePerResidueScores(unordered_map<Residue*, mstreal> scores);

/**
 Virtual class that scores each residue in a seed structure.
 */
class SeedScorer {
public:
    SeedScorer(Structure *target): _target(target) {}
    
    /**
     Score each residue in a seed structure.
     
     @param seed the seed Structure, not including chains that align with the
            target
     */
    virtual unordered_map<Residue *, mstreal> score(Structure *seed) = 0;
    
    
protected:
    Structure *_target;
};

/**
 Virtual SeedScorer subclass that provides useful functionality for interfacing
 with a FASST database.
 */
class FASSTScorer: public SeedScorer {
public:
    /**
     Initializes a FASSTScorer and loads the FASST database.
     
     @param target the target structure, whose sequence is known
     @param configFilePath the path to the binary FASST database
     @param fractionIdentity the redudancy cutoff when searching for matches
     @param maxNumMatches the number of matches after which to stop searching
     @param vdwRadius NOT the actual VDW radius within which clashes are checked,
            but the *ratio* of the maximum VDW radius for the specific atom within
            which clashes are checked
     @param scoreType the score type
     */
    FASSTScorer(Structure *target, string configFilePath, double fractionIdentity = 0.4, int maxNumMatches = 8000, double vdwRadius = 1.0);
    FASSTScorer(Structure *target, string configFilePath, int scoreType, double fractionIdentity = 0.4, int maxNumMatches = 8000, double vdwRadius = 1.0);
    ~FASSTScorer();
    
    virtual unordered_map<Residue *, mstreal> score(Structure *seed) = 0;
    
    /**
     Checks for steric clashes between the seed and the target, and returns false
     if they clash. If not, cleans up the seed structure, and appends the
     seed chain into targetStructBB.
     
     @param seed the seed structure
     @return true if the combined structure was created, and false if there was
     a steric clash
     */
    bool prepareCombinedStructure(Structure *seed);
    bool prepareCombinedStructure(Structure *seed, bool clearNames);

    /**
     Removes the seed chain from targetStructBB.
     */
    void resetCombinedStructure();
    
    /*
     Toggle whether the scorer allows clashes
     */
    void setIgnoreClash(bool _ignoreClash) {ignoreClash = _ignoreClash;}
    
    void setQueryHomologyCutoff(double val) {
        if (val <= 0.0) MstUtils::error("Query homology cutoff cannot be <= 0","FASSTScorer::setQueryHomologyCutoff");
        else queryHomologyCutoff = val;}
    
protected:
    configFile config;
  
    FASST* fasst;
    RotamerLibrary *rl;
    
    double fractionIdentity;
    int maxNumMatches;
    double vdwRadius;
    bool ignoreClash = false;
    
    /**
     Any match with a sequence identity greater than this cutoff will be removed from the matches, as it is considered homologous to the
     query and thus redundant. Any value above 1.0 signals that this should not be applied when searching
     */
    double queryHomologyCutoff = 2.0;
    
    // Backbone-only target
    Structure targetStructBB;
    set<Residue*> targetResidueSet;
    
    // Proximity search
    AtomPointerVector psTargetAPV;
    ProximitySearch ps;
    
    void loadFromTarget(string fasstDB);
    void loadFromTarget(string fasstDB, int scoreType);
    
    /**
     Remaps the residue keys from the given temporary result so that they are
     represented by residues at the appropriate indexes in seed.
     
     @param tempResult a score map with keys corresponding to targetStructBB
     @param seed a structure containing residues to remap to
     @return a remapped dictionary of scores
     */
    unordered_map<Residue*, mstreal> remapResiduesFromCombinedStructure(unordered_map<Residue*, mstreal> tempResult, Structure *seed);
    
    /**
     Convenience method that returns a map of each residue in seed to DBL_MAX.
     */
    unordered_map<Residue *, mstreal> invalidScoreMap(Structure *seed);
};

/*
 A helper class that keeps track of the number of seeds that each target residue contacts.
 
 The class is designed to be used in parallel when calculating the total number of contacts, as this
 can be relatively slow. The seeds can be divided into batches, which are assigned to separate instances
 of contactCounter. When all are complete, the contacts are written to contactFiles, and then are read
 in such that the results are combined and can be rewritten as a single file. See programs/countContacts.cpp
 */

class contactCounter {
public:
    /**
     Constructs the counter with zero counts for each target residue structure
     
     @param _target structure of the target protein (can be more than one chain)
     @param _RL pointer to a rotamer library
     @param _cParams the contact parameters
     @param _flankingRes the number of residue on the N/C terminus of the seeds to ignore
     */
    contactCounter(Structure* _target, RotamerLibrary *_RL, contactParams _cParams, int _flankingRes = 2) : target(_target), RL(_RL), cParams(_cParams), flankingRes(_flankingRes) {
        buildContactCounts();
        targetCopy = new Structure(*target);
    };
    
    /**
     Constructs the counter with counts for each target residue structure specified in the contactsFile
     
     @param contactsFile path to file with number of seed contacts to each residue in the target
     @param _target structure of the target protein (can be more than one chain)
     @param _RL pointer to a rotamer library
     @param _cParams the contact parameters
     @param _flankingRes the number of residue on the N/C terminus of the seeds to ignore
     */
    contactCounter(string contactsFile, Structure* _target, RotamerLibrary *_RL, contactParams _cParams, int _flankingRes = 2) : target(_target), RL(_RL), cParams(_cParams), flankingRes(_flankingRes) {
        buildContactCounts();
        readContactsFile(contactsFile);
        targetCopy = new Structure(*target);
    };
    
    ~contactCounter() {
        delete targetCopy;
        if (contact_out != nullptr) delete contact_out;
    };
    
    void setWriteSeedContacts(string contactsFile) {
        contact_out = new fstream();
        MstUtils::openFile(*contact_out,contactsFile,fstream::out);
        *contact_out << "seed_name,chain_id,residue_num,residue_idx,num_contacts" << endl;
    }
    
    /**
     Counts the number of unique contacts between the target protein and the provided seed and
     adds this to contactCounts
     
     @param seed the structure of the seed (must be a single chain and have a complete backbone)
     @return the number of new contacts between the target and seed
     */
    int countContacts(Structure *seed);
    
    /**
     Gets the number of contacts to a specific residue in the target protein
     
     @param res a residue in target
     @return the number of contacts between res and observed seeds
     */
    int getContactCounts(Residue *res) {
        if (contactCounts.count(res) != 1) MstUtils::error("Provided residue does not exist in target structure","contactCounter::getContactCounts");
        return contactCounts[res];
    }
    
    /**
     Sets the contact count for each residue to 0
     */
    void reset() {
        for (auto it : contactCounts) it.second = 0;
    };
    
    /**
     Reads the contacts from a contactsFile and adds these to current contact counts
     
     @param contactsFile path to file with number of seed contacts to each residue in the target
     */
    void readContactsFile(string contactsFile);
    
    /**
     Writes the current contacts count to contactsFile
     
     @param contactsFile path to file with number of seed contacts to each residue in the target
     */
    void writeContactsFile(string contactsFile);
    
protected:
    void buildContactCounts() {
        for (Residue* R : target->getResidues()) contactCounts[R] = 0;
    }
    
//    vector<Residue*> getSeedRes(Chain* seedChain) {
//        vector<Residue*> seedRes;
//        if (seedChain->residueSize() > flankingRes*2) {
//            for (int i = flankingRes; i < seedChain->residueSize() - flankingRes; i++) seedRes.push_back(&seedChain->getResidue(i));
//        }
//        return seedRes;
//    }
    
    Residue* getTargetResidue(Residue* targetCopyRes) {
        /*
         since the seed chain is always added as the final chain, the index of the residues
         in the previous chains (i.e. the target) should not change
         */
        return &target->getResidue(targetCopyRes->getResidueIndex());
    }
    
    vector<Residue*> eliminateFlankingRes(vector<Residue*> allRes) {
        vector<Residue*> filtRes;
        for (Residue* R : allRes) {
            int idxInChain = R->getResidueIndexInChain();
            if (min(int(allRes.size())-1-idxInChain,idxInChain-0) >= flankingRes) {
                // if source residue is not less than flankRes from either termini, add
                filtRes.push_back(R);
            }
        }
        return filtRes;
    }
    
private:
    Structure *target;
    Structure *targetCopy; //copied structure is modified when counting contacts
    map<Residue*, int> contactCounts; //Residue pointer from target
    int flankingRes;
    contactParams cParams;
    RotamerLibrary *RL;
    bool verbose = true;
    fstream* contact_out = nullptr;
};

/**
 Scores each residue in a seed structure according to the sequence-structure
 compatibility score, defined in Craig Mackenzie's thesis, section 4.3.2. Adapted
 from structgen, score.h/.cpp.
 
 To use the new formulation of sequence score, which includes a "contact score"
 to normalize for the probability of observing a particular seed structure given
 a target structure, the method collectContacts must be called with every seed
 before starting to score. These results can be written to file and read back
 with the writeContactCounts and readContactCounts methods.
 
 */
class SequenceStructureCompatibilityScorer: public FASSTScorer {
public:
    
    /**
     Initializes a sequence compatibility scorer and loads the FASST database.
     
     @param target the target structure, whose sequence is known
     @param rParams the parameters to use for determining the RMSD cutoff for
            search
     @param contParams the criteria for determining if a contact exists
     @param configFilePath the path to the FASST database to load
     @param dTERMenConfigFilePath the path to the dTERMen config file to load
     @param targetFlank the number of flanking residues on either side of the
            target residue
     @param seedFlank the number of flanking residues on either side of the seed
            residue
     @param fractionIdentity the redudancy cutoff when searching for matches with
            FASST
     @param minRatio the minimum ratio of pair-TERM probability to self-TERM
            probability; if smaller, the residue's score will be set to inf
     @param pseudocount the pseudocount to use when generating sequence statistics
            for each position in the fragments
     @param minNumMatches the minimum number of matches for the pair-TERM needed
            to compute the score; otherwise the residue's score is set to inf
     @param maxNumMatches the maximum number of matches in FASST before the
            search terminates
     @param vdwRadius NOT the actual VDW radius within which clashes are checked,
            but the *ratio* of the maximum VDW radius for the specific atom within
            which clashes are checked
     */
    SequenceStructureCompatibilityScorer(Structure *target, rmsdParams& rParams, contactParams& contParams, string configFilePath, string dTERMenConfigFilePath, bool filter = false, mstreal freedomLimInput = 0.5, int scoreTypeInput = 1, int targetFlank = 2, int seedFlank = 2, double fractionIdentity = 0.4, double minRatio = 0.05, double pseudocount = 0.25, int minNumMatches = 1, int maxNumMatches = 8000, double vdwRadius = 1.0);
    ~SequenceStructureCompatibilityScorer();
    
    /**
     This currently only supports seeds with a single chain.
     
     @return a map containing all residues for which a contact was found, and
            their scores. If the seed clashes with the target, or no contacts
            were found, an empty map is returned.
     */
    unordered_map<Residue *, mstreal> score(Structure *seed) override;
    
    /*
     Collects contacts from the given seed structure without scoring them. This
     must be done for all seeds before scoring them, or else the contact score
     cannot be computed.

     @param seed the seed structure from which to compute and count contacts
     */
    void collectContacts(Structure *seed);

    /**
     Writes the contact counts for each target residue to the given file path
     as CSV.

     @param filePath the path to write to
     */
    void writeContactCounts(string filePath);

     /*
     Reads the contact counts for each target residue from the given CSV file.
     
     @param filePath the path to read from
     */
    void readContactCounts(string filePath);
    
    /**
     Calculates the fraction of total contacts with the target residue that
     overlap with the seed residue. Requires that (1) this scorer has access to
     an adjacency graph, and (2) all residues that could potentially contact
     targetRes have passed through a score() call. This method works even if
     noQueries is false, however.
     
     @param seedRes the seed residue with which to search for overlaps
     @param targetRes the residue on the target within which to find contacts
     @param pseudocount a value added to both the number of overlaps and the
            number of contacts before finding the ratio between them
     @return the contact score, i.e. the probability of the given seed residue
            position conditioned on the target residue position
     */
    mstreal contactScore(Residue *seedRes, Residue *targetRes, double pseudocount = 1);
    
    /**
     The workhorse of this class, this method computes the sequence-structure
     score as adapted from structgen.
     
     @param seed the original seed structure
     @param combStruct a structure containing both the target and the seed
     @param cl a list of contacts to check
     @param toScore the residues on the target that are being scored
     @return a map containing all residues for which a contact was found, and
            their scores. Any search in which the number of results is
            insufficient or the probability ratio is too small results in the
            seed residue's score being set to infinity.
     */
    unordered_map<Residue*, mstreal> sequenceStructureScore(Structure *seed, Structure& combStruct, contactList& cl, set<Residue*>& toScore, ConFind &confindCmbined, ConFind &confindTarget);

    /**
     Computes the sequence-structure score of a single contact. This method is
     pulled out from the main sequenceStructureScore method so that a different
     scoring method can easily apply the same subroutines used in this scorer.
     
     @param frag the fragment containing the two contacting residues and their
            flanking regions (expansions)
     @param scoringRes the residue on the target to be scored
     @param seedRes the contacting residue on the seed
     @param seenProbs a cache of fragment probabilities that, if non-null, can
            be used to skip the computation
     @param cacheBackground whether to read/write cached values for the
            background probabilities (probability of self-TERM around scoringRes)
    @param confindCombind relating scoring residues to seed residues
     @return the sequence-structure score of the given contact, or DBL_MAX
            (infinity) if there were not enough searfch results, or if the
            probability ratio is below minRatio
     */
    mstreal sequenceStructureScoreComponent(Fragment& frag, Residue* scoringRes, Residue* seedRes, ConFind &confindCombined, map<Fragment, double> *seenProbs = nullptr, bool cacheBackground = true);

    /// Path at which to write out queries.
    string *queryWritePath = nullptr;
    
    /// Path at which to write out score breakdowns for each contact.
    string *scoresWritePath = nullptr;

    /// Path at which to write out which residues were filtered.
    string *filterWritePath = nullptr;

    /// Path at which to write out freedom info for matches.
    string *freedomWritePath = nullptr;
        
    /// Stores the scores for each residue and its neighborhood.
    SeedGraphMap<mstreal> *adjacencyGraph = nullptr;

    /// If true, then don't actually calculate the scores, just log that they were "calculated"
    bool noQueries = false;

    // Max score
    mstreal maxScore = 100;
    
    // Analytics
    
    /// @return the number of total (seed) residues that this scorer has scored, de novo or not
    int numResiduesScored() { return residuesScored; }
    /// @return the number of unique residues that this scorer has scored
    int numUniqueResiduesScored() { return uniqueResiduesScored; }
    /// @return the number of total seed queries that have been made (i.e. unique residues scored)
    int numSeedQueries() { return seedQueries; }
    /// @return the number of total target residue queries that have been made
    int numTargetQueries() { return targetQueries; }
    /// @return the number of matches found for each query
    vector<mstreal> numMatches() { return matches; }


    
private:
    FragmentParams fragParams;
    rmsdParams rParams;
    contactParams contParams;
    string dTERMenConfigFile;
    int scoreType;
    double minRatio;
    double pseudocount;
    double vdwRadius;
    int minNumMatches;
    bool mustContact = true;
    int targetFlank;
    int seedFlank;
    
    // If true, filter residues to score based on freedom
    bool freedomFilter = false;
    double freedomLim;

    /// If true, the seeds will not be scored; instead the contacts with target will be counted
    bool countingContacts = false;
    
    // Analytics
    int residuesScored = 0;
    int uniqueResiduesScored = 0;
    int seedQueries = 0;
    int targetQueries = 0;
    vector<mstreal> matches;

    Structure *targetOnlyStructure;
    
    /**
     Cache the background probabilities for each aa in target (pair includes
     total # residues in flanking region for verification, and probability).
     */
    unordered_map<Residue *, pair<int, mstreal>> backgroundProbs;
    
    /// Number of contacts analyzed for each residue in the target.
    unordered_map<Residue *, int> targetContactCounts;
    
    /// Write stream for writing out scores
    ofstream *scoreWriteOut = nullptr;

    // Write stream for writing out filtered residues
    ofstream *filterWriteOut = nullptr;

    // Write stream for writing out freedoms
    ofstream *freedomWriteOut = nullptr;
    
    /**
     Initializes the score write stream if a path was provided, and returns it.
     Returns null if scoresWritePath is null.
     */
    ofstream *getScoreWriteStream();

    /**
     Initializes the filter write stream if a path was provided, and returns it.
     Returns null if filterWritePath is null.
     */
    ofstream *getFilterWriteStream();

    /**
     Initializes the freedom write stream if a path was provided, and returns it.
     Returns null if freedomWritePath is null.
     */
    ofstream *getFreedomWriteStream();
    
    /**
     Trims the given fragment around the two central residues so that the
     targetFlank and seedFlank parameters are obeyed. Places the resulting
     structure in result.
     */
    vector<Residue *> trimFragment(Fragment &frag, Residue *scoringRes, Residue *seedRes, Structure &result);
};

/**
 Scores a seed structure using the structural compatibility score, defined in
 Craig Mackenzie's thesis, section 4.3.1. Adapted from structgen, score.h/.cpp.

 This score is essentially the number of designable contacts per residue in the
 seed chain - in this implementation, higher is better. If a residue has a non-
 designable fragment, its score is set to -infinity.
 
 note: if scoreAll is set to true, then non-designable contacts are not added
 to the final score.
 
 renamed from StructureCompatibilityScorer to SeedDesignabilityScorer
 */
class SeedDesignabilityScorer: public FASSTScorer {
public:
    /**
     For descriptions of parameters, see initializer for
     SequenceCompatibilityScorer.
     
     option added by sebastian on 20/07/19
     @param scoreAll search all contacts, regardless if some are found to be non-designable
     */
    SeedDesignabilityScorer(Structure *target, FragmentParams& fragParams, rmsdParams& rParams, contactParams& contParams, string configFilePath, double fractionIdentity = 0.4, int minNumMatches = 1, int maxNumMatches = 8000, double vdwRadius = 0.7, bool _scoreAll = false): FASSTScorer(target, configFilePath, fractionIdentity, maxNumMatches, vdwRadius), fragParams(fragParams), rParams(rParams), contParams(contParams), minNumMatches(minNumMatches), scoreAll(_scoreAll) {};
    ~SeedDesignabilityScorer() {
    }
    
    /**
     This currently only supports seeds with a single chain.
     
     @return a map containing all residues for which a contact was found, and
             their scores. If the seed clashes with the target, or no contacts
             were found, an empty map is returned.
     */
    unordered_map<Residue *, mstreal> score(Structure *seed) override;
    
    /*
     Option added by sebastian 20/07/19
     @param intra when true, defines contacts between residues of the path
     */
    
    void score(Structure *seed, mstreal &totalScore, int &numContacts, int &numDesignable, bool intra = false, bool score_all = false);

    bool clashes(Structure *seed);
    
    void setScoreAll(bool _scoreAll) { scoreAll = _scoreAll;}

    /**
     The workhorse of this class - this method counts the number of designable
     contacts in combStruct for each residue in seedResidues.
     
     @param combStruct the target and seed structures together
     @param cl a list of contacts between target and seed residues
     @param seedResidues a list of seed residues, used to reference which member
            of a contact should be scored
     @return a map containing the number of designable contacts for each residue
             in seedResidues
     */
    unordered_map<Residue *, mstreal> designabilityScore(Structure& combStruct, contactList& cl, vector<Residue*> seedResidues, string seedName = "");
    
//    /// Path at which to write out queries.
//    string *queryWritePath = nullptr;
    
    /// Path at which to write out score breakdowns for each contact.
    string scoresWritePath = "";
    
private:
    FragmentParams fragParams;
    rmsdParams rParams;
    contactParams contParams;
    double vdwRadius;
    int minNumMatches;
    bool mustContact = true;
    bool scoreAll;

    // Statistics for current seed
    int _numDesignable = 0;
    
    /// Write stream for writing out scores
    ofstream *scoreWriteOut = nullptr;
    
    /**
     Initializes the score write stream if a path was provided, and returns it.
     Returns null if scoresWritePath is null.
     */
    ofstream *getScoreWriteStream();
};



#endif /* seedscore_h */
