#ifndef _DTERMENSCORER_H
#define _DTERMENSCORER_H

#include "msttypes.h"
#include "mstfasst.h"
#include "mstsequence.h"
#include "mstcondeg.h"
#include "mstmagic.h"
#include "dtermen.h"
using namespace MST;

class EnergyTable;
class termData;


class dTERMenScorer {
  public:
    dTERMenScorer();
    dTERMenScorer(const string& configFile);
    void init();
    void readConfigFile(const string& configFile);
    void setEnergyFunction(const string& ver); // sets the energy function version and alters any necessary parameters
    void setRecordFlag(bool record = true) { recordData = record; }

    mstreal getkT() const { return kT; }
    FASST* getFASST() { return &F; }
    RotamerLibrary* getRotamerLibrary() { return &RL; }
    vector<res_t> getGlobalAlphabet() { return globAlph; }
    void setAminoAcidMap();
    void printAminoAcidMap();
    int globalAlphabetSize() const { return globAlph.size(); }
    mstreal getHomologyCutoff() const { return homCut; }

    bool isInGlobalAlphabet(const string& aa) const;
    bool isInGlobalAlphabet(res_t aa) const;
    // the following four functions convert to and from the "internal" amino-acid index
    int aaToIndex(const string& aa) const;
    int aaToIndex(res_t aa) const;
    bool aaIndexKnown(int aaIdx) const;
    string indexToResName(int idx) const;
    res_t indexToAA(int idx) const;

    struct oneDimHist {
      vector<vector<int> > bins;
      vector<mstreal> binEdges;
      vector<mstreal> binMasses;
      vector<mstreal> weights;
    };
    struct twoDimHist {
      vector<vector<vector<int> > > bins;
      vector<mstreal> xBinEdges, yBinEdges;
      vector<mstreal> weights;
    };
    struct zeroDimPotType { // NOTE: in the future, these should become classes with overloaded "access" operators for lookup
      vector<mstreal> aaEnergies;
    };
    struct oneDimPotType {
      vector<mstreal> binEdges;
      vector<vector<mstreal> > aaEnergies;
    };
    struct twoDimPotType {
      vector<mstreal> xBinEdges;
      vector<mstreal> yBinEdges;
      vector<vector<vector<mstreal> > > aaEnergies;
    };
    void buildBackgroundPotentials();
    zeroDimPotType buildZeroDimPotential(const vector<int>& AA, vector<vector<mstreal> >& backPot);
    oneDimPotType buildOneDimPotential(const oneDimHist& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot = false);
    oneDimPotType buildOneDimPotential(const oneDimHist& H, const vector<int>& AA, mstreal pc); // without specifying a background potential
    twoDimPotType buildTwoDimPotential(const twoDimHist& H, const vector<int>& AA, mstreal pc, vector<vector<mstreal> >& backPot, bool updateBackPot = false);
    mstreal lookupZeroDimPotential(const zeroDimPotType& P, int aa);
    mstreal lookupOneDimPotential(const oneDimPotType& P, mstreal x, int aa);
    mstreal lookupTwoDimPotential(const twoDimPotType& P, mstreal x, mstreal y, int aa);
    void printZeroDimPotential(const zeroDimPotType& P);
    void printOneDimPotential(const oneDimPotType& P);
    void printTwoDimPotential(const twoDimPotType& P);

    /* Bins the data in the input vector X according to the binning type and
     * parameters. Outputs a struct with members bins and binEdges. The k-th bin
     * (k being between 0 and n-1, where n is the number of bins) is defined as
     * the interval [binEdges[k]; binEdges[k+1]), except for the last bin, which
     * does include the right limit: [binEdges[k]; binEdges[k+1]]. Upon returning,
     * bins[k] will contain the list of indices of data points that map into the
     * k-th bin. Supported binning types are (values for binSpecType):
     * 1 -- uniform binning. binSpec is expected to be {min value, max value,
     *      and desired number of bins}.
     * 2 -- non-uniform binning with some minimal number of elements per bin and
     *      a minimal bin width. binSpec is expected to be {min value, max value,
     *      minimum number of points per bin, minimal bin width}.
     * Optional argument M (must be the same size as X, if specified) supplies
     * the multiplicity of each point. This is used in non-uniform binning to
     * count data "mass" (i.e., sum of inverses of multiplicities) rather the
     * pure number of points. If isAngle is set to true, will treat input data
     * as angles in degrees, mapping them to the interval [-pi; pi). This makes
     * it so that the bin definitions above do the right thing of counting +/- pi
     * just once.
     * NOTE: any data points that fall outside of the range [min; max] are
     * ignored in building the potential. This may include bad diehdral angles,
     * such as phi/psi angles for terminal residues. */
    oneDimHist binData(const vector<mstreal>& X, int binSpecType, const vector<mstreal>& binSpec, const vector<mstreal>& M = vector<mstreal>(), bool isAngle = false);
    twoDimHist binData(const vector<mstreal>& X, const vector<mstreal>& Y, const vector<mstreal>& xBinSpec, const vector<mstreal>& yBinSpec, const vector<mstreal>& M = vector<mstreal>(), bool isAngle = false);

    void readBackgroundPotentials(const string& file);  // TODO
    void writeBackgroundPotentials(const string& file); // TODO
    void writeRecordedData(const string& file); // writes all recorded TERM data to a file

    mstreal backEner(const string& aa) { return lookupZeroDimPotential(bkPot, aaToIndex(aa)); }
    mstreal bbOmegaEner(mstreal omg, const string& aa) { return lookupOneDimPotential(omPot, omg, aaToIndex(aa)); }
    mstreal bbPhiPsiEner(mstreal phi, mstreal psi, const string& aa) { return lookupTwoDimPotential(ppPot, phi, psi, aaToIndex(aa)); }
    mstreal envEner(mstreal env, const string& aa) { return lookupOneDimPotential(envPot, env, aaToIndex(aa)); }

    /* Compute the dTERMen self energy of a given Residue; identifies relevant
     * contacts using ConFind. The first form computes the value for a specific
     * amino acid at the given position (by default the one actually there), the
     * second form returns values for all amino acids at the position, and the
     * third form works as the second, but accepts a pre-built ConFind object.
     * When multiple energies are needed from the same structure and, passing the
     * same ConFind object will make it so that contacts and freedoms will only
     * need to be computed once. */
    mstreal selfEnergy(Residue* R, const string& aa = "");
    vector<mstreal> selfEnergies(Residue* R, bool verbose = false);
    vector<mstreal> selfEnergies(Residue* R, ConFind& C, bool verbose = false, bool seedIndicator = false);

  protected:
    int findBin(const vector<mstreal>& binEdges, mstreal x); // do a binary search to find the bin into which the value falls
    mstreal backEner(int aai) { return lookupZeroDimPotential(bkPot, aai); }
    mstreal bbOmegaEner(mstreal omg, int aai) { return lookupOneDimPotential(omPot, omg, aai); }
    mstreal bbPhiPsiEner(mstreal phi, mstreal psi, int aai) { return lookupTwoDimPotential(ppPot, phi, psi, aai); }
    mstreal envEner(mstreal env, int aai) { return lookupOneDimPotential(envPot, env, aai); }
    void printSelfComponent(const CartesianPoint& ener, const string& prefix);

    /* Given a list of FASST solutions, computes the residual statistical energy
     * for all amino acids at the position with index cInd, after accounting for
     * all "trivial" background statistical contributions at this position
     * accross all of the matches. */
    CartesianPoint singleBodyStatEnergy(fasstSolutionSet& matches, int cInd, int pc);

    /* For a given match, compute the expectation of any given amino acid at the
     * specified position based on "trivial" background statistical contributions. */
    CartesianPoint backExpectation(const fasstSolution& m, int cInd);

    /* Counts observations at the given positions across all matches. */
    CartesianPoint singleBodyObservations(fasstSolutionSet& matches, int cInd);

    /* Same as above, but for amino-acid observations at two sites. */
    CartesianPoint twoBodyObservations(fasstSolutionSet& matches, int cIndI, int cIndJ);

    /* Computes the number of times each amino acid is expected to be found at
     * the given position, across all matches. NOTE: skips any matches that have
     * an unknown residue in the position. This is because when comparing to the
     * observed counts, such matches would have no "vote", so it would be unfair
     * to given them a vote in computing the expectation. The optional pointer
     * can be specified to collect the expectation in each match. */
    CartesianPoint singleBodyExpectations(fasstSolutionSet& matches, int cInd, vector<CartesianPoint>* breakDown = NULL);

    /* Computes the number of times each amino acid is expected to be found at
     * the given position in each match. Identifies an underlying amino-acid bias
     * energy vector to make sure that the marginals (i.e., the total number of
     * times each amino acid is expected across all matches) approximately equal
     * to the expectations. In practice, the agreement should be essentially
     * perfect, limited only by the number of iterations of the underlying itera-
     * tive procedure. */
    vector<CartesianPoint> singleBodyExpectationsMatchedMarginals(fasstSolutionSet& matches, int cInd);

    mstreal enerToProb(vector<mstreal>& ener);
    mstreal enerToProb(const vector<mstreal>& _ener) { vector<mstreal> ener = _ener; return enerToProb(ener); }

    /* Amino-acid indices provide a convenient way to index into a array, while
     * also encoding the amino-acid type. For amino-acid pairs, this is tricky,
     * because a 2D array is needed. These functions convert between two amino-
     * acid indices and a single index that uniquely identifies the pair. */
    int pairToIdx(int aai, int aaj) const;
    pair<int, int> idxToPair(int idx) const;

    /* Given a list of source residues, and a ConFind object (already initialized
     * with the relevant structure), returns a list of residue pairs that corres-
     * pond to contacts with source residues. These are either intra contacts--,
     * contacts formed between the source residues--or contacts between source
     * residues and residues not belonging to that list. This is controlled by the
     * last integer parameter (0, the default, corresponds to inter, 1 to intra,
     * and two to both). In the case of inter-contacts, the source residue of
     * each contact will be listed first in each pair (i.e., will be stored in
     * the "first" field of each pair). The contacts in the returned list will
     * be ordered as follows: first contact-degree based contacts, ordered by
     * contact degree in descending order, then any additional contacts emerging
     * from interference, ordered by interference. */
    vector<pair<Residue*, Residue*>> getContactsWith(const vector<Residue*>& source, ConFind& C, int type = 0, bool verbose = false);

  private:
    FASST F;
    fasstSearchOptions foptsBase; // base FASST options that will be used with every search
    RotamerLibrary RL;
    string fasstdbPath, backPotFile, rotLibFile, efunVer;
    zeroDimPotType bkPot;
    oneDimPotType omPot, envPot;
    twoDimPotType ppPot;
    mstreal kT, cdCut, intCut, selfResidualPC, selfCorrPC, homCut;
    int pmSelf, pmPair;
    int selfResidualMinN, selfResidualMaxN, selfCorrMinN, selfCorrMaxN, selfCorrMaxCliqueSize, pairMinN, pairMaxN;
    bool recordData;
    vector<termData> data;
    Sequence targetOrigSeq;
    vector<int> variableResidues;
    map<string, vector<mstreal>> targetResidueProperties;

    /* We may want to deal with different "universal" alphabets (separate from
     * the design alphabet). We may want to interpret SEC (selenocysteine), for
     * example, as CYS (cysteine) in gathering sequence statistics. Or, we may
     * want to keep these as separate counts. The following several variables
     * support this capability by providing a mapping between all the possible
     * amino acid names and indices defined in SeqTools to a smaller sub-set,
     * defined over a contiguous set of indices starting with 0. In all cases
     * residues are identified by their res_t index, and the classification
     * between residue strings and res_t is left to SeqTools. But, for simplicity,
     * in the comments below, I will refer to residues by their string.
     * aaMap -- stores any mapping between any non-standard amino acid names and
     *          standard ones. E.g., aaMap["SEC"] may contain "CYS", saying that
     *          SEC is interpreted as a CYS, or it may not exist, which says that
     *          we want to explicitly track the counts of SEC, separate from CYS.
     * globAlph -- a list of all counted amino acids, referred to by their most
     *             "standard" name. E.g., if aaMap["SEC"] == "CYS", then globAlph
     *             will contain "CYS" but not "SEC". Further, amino-acid names
     *             in globAlph are stored according to their contiguous index
     *             (internal ones for dTERMen, not the same as SeqTools indices).
     * aaIdx    -- effectively the opposite of globAlph. For every counted amino
     *             acid, aaIdx will contain its corresponding internal index. E.g.,
     *             aaIdx["SEC"] and aaIdx["CYS"] would have the same index, if
     *             if aaMap["SEC"] == "CYS".
     * aaMapType -- specifies the type of mapping between non-standard amino
     *              acids and their standard counterparts. Used by setAminoAcidMap() */
    map<res_t, res_t> aaMap;
    vector<res_t> globAlph;
    map<res_t, int> aaIdx;
    int aaMapType;

};

#endif