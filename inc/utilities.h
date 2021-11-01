#ifndef utilities_h
#define utilities_h

#include <list>
#include <regex>
#include <unordered_map>

#include "msttypes.h"
#include "mstcondeg.h"
#include "mstlinalg.h"
#include "mstsequence.h"

#include "mstsystem_exts.h"
#include "vdwRadii.h"

using namespace MST;

/* A simple class that stores paths to a rotamer library/fasst database.
 
 These can vary between systems, so I found it easier to maintain a configFile for each project
 and call that when needed.
 */
class configFile {
public:
    configFile(string configFile_path);
    
    string getRL() {return RL_path;}
    string getDB() {return fasst_DB_path;}
protected:
private:
    string RL_path;
    string fasst_DB_path;
};


class generalUtilities {
    
public:
    static string selectionStringFromContactingRes(contactList conts);
    static string selectionStringFromRes(vector<Residue*> res);
    
    // copy the coordinates of one structure onto another
    static void copyAtomCoordinates(Structure* S_dest, const Structure* S_source);
    static void copyAtomCoordinates(Structure* S_dest, vector<Atom*> atoms_source);
    
    // Import structures from a binary file
    static vector<Structure*> readStructuresBin(string binFile);
    static vector<Structure*> readStructuresList(string listFile);
    static void writeStructuresBin(vector<Structure*> structures, string binFile);
    static void writeStructuresBin(vector<Structure> structures, string binFile);
    
    static vector<pair<Residue*, Residue*>> getContactsWith(const vector<Residue*>& source, ConFind& C, int type, mstreal cd_threshold, mstreal int_threshold, mstreal bbInteraction_cutoff, bool verbose);
    static vector<Residue*> getContactingResidues(vector<pair<Residue*,Residue*>> cont_pairs);
    
    static contactList getContactsWith(const vector<Residue*>& source, ConFind& C, mstreal threshold, string type);
    static set<pair<Residue*,Residue*>> mergeContactLists(const vector<contactList> CL);
    
    // T should be a vector containing some objects
    template <class T>
    static T generateKCombinations(T initial_set, int k);
    
    /* recursive algorithm for generating all unique combinations of k residues from a larger set (without replacement).
     
     initial set - the set of residues that the user wants to generate combinations from
     k           - the number of elements in the combinations
     
     The user should envision a tree structure where each node represents an initial set with some
     associated k. For each member of the set the function will be called recursively with 1) k
     decremented by 1 and 2) the selected member removed from the set. For each recursive call to the
     function, a new node, one level below is generated with the associated variables and an edge is drawn.
     At the bottom of the tree, an empty vector of vectors is returned. If the initial set becomes empty
     before the bottom of the tree is reached, the entire branch becomes void.
     
     Example:
     
     from [1,2,3] choose all combinations of two objects
     
     k = 2               [1,2,3]
     current res       1/  2|  3\
     k = 1          [2,3]  [3]  []
     current res   2/ 3\   3|
     k = 0       [3]   []  []
     
     each branch that reaches the level where k = 0 generates a combination; the current residues
     along each path from the top node to the bottom nodes form the set.
     */
    
    static vector<vector<Residue*>> generateAllCombinationsKRes(vector<Residue*> initial_set, int k);
    
    /* A modified version of the above function that returns combinations for all values of k,
     between 1 and the size of the input vector.
     
     This is done by simple calling the above function, for each k (in descending order).
     
     */
    
    static vector<vector<Residue*>> generateAllCombinationsRes(vector<Residue*> initial_set);
    
    
    // The following functions were adapted from msttypes (originally developed by Craig for the set cover library)
    // (L0 was determined to be around 15)
    /* L stores the number of residues in each chain, assumed to be contiguous.
     * Residues in different chains are assumed to be statistically independent. */
    static mstreal fragDegreesOfFreedom(const vector<int>& L, mstreal L0 = 15);
    
    /* I stores residue indices of residues in each chain, so that non-contiguous
     * residues within a chain can be treated. Residues in different chains are
     * still assumed to be independent. */
    static mstreal fragDegreesOfFreedom(const vector<vector<int> >& I, mstreal L0 = 15);
    
    /* Takes every chain in S to be independent. */
    static mstreal fragDegreesOfFreedom(const Structure& S, mstreal L0 = 15);
    
    /* DOF is computed for the fragment comprising residues from the given
     * Structure S, whose indices are given in J. Like in the previous function,
     * for residues located within the same chain of S, their relative locations
     * are accounted for when computing degrees of freedom. */
    static mstreal fragDegreesOfFreedom(const vector<int>& J, const Structure& S, mstreal L0 = 15);
    
    /*
     Defines a vector normal to the plane formed by N, Ca, and C atoms in each residue and
     computes the cosine of the angle between them. Since this function just compares two
     vectors defined in the same way, I don't think it matters how they are defined precisely
     */
    static mstreal cosAngleBetweenNormalVectors(Residue* R1, Residue* R2);
    
    static mstreal cosAngle(CartesianPoint v1, CartesianPoint v2);
    
    static mstreal avgCosAngleBetweenSegments(const vector<Residue*> &seg1, const vector<Residue*> &seg2);
    
    /*
     Verifies that the residues are covalently bonded in the order that they are provided
     */
    static bool contiguousResidues(vector<Residue*> segment, mstreal maxPeptideBond = 2.0);
    
    // Uses brute-force to find the best alignment (no gaps) of two peptide structures
    static mstreal bestRMSD(Chain* C1, Chain* C2);

};

class volumeCalculator {
public:
    volumeCalculator(Structure* _protein, mstreal _qR = 12.0) : protein(_protein), qR(_qR) {
        atoms = protein->getAtoms();
        ps = ProximitySearch(atoms,qR/2);
    };
  
    volumeCalculator() {};
    
    void setStructure(Structure* _protein) {
        protein = _protein;
        atoms = protein->getAtoms();
        ps = ProximitySearch(atoms,qR/2);
    }

    void setQueryRadius(mstreal _qR) {
        if (atoms.size() == 0) MstUtils::error("Structure must be set","volumeCalculator::setQueryRadius");
        qR = _qR;
        ps = ProximitySearch(atoms,qR/2);
    }

    /**
     An analytical approach to calculating the fraction of volume around a residue occupied by atoms.
     
     NOTE: does not work. In order to accurately calculate the occupied volume using this approach must consider pairwise and higher
     order overlaps, otherwise will double count. This is hard to do, see this paper that used complicated geometry to solve this:
     https://pubs.acs.org/doi/pdf/10.1021/ja00291a006
     
     The idea was to sum up the volume of VDW spheres within the "query sphere", the space we were searching for atoms. From this we
     would subtract the pairwise overlaps between VDW spheres, which otherwise would be double counted. The issue is that this
     overcorrects, so in cases where we have higher-order overlaps (which does happen) we would be missing that volume. Beyond this,
     when finding pairwise overlaps, it is important to calculate how much of this falls inside the query sphere, but this is also difficult and
     not done here. Instead that area is just ignored.
     
     @param R the residue around which we will do the calculation
     @return the fraction of the sphere around the alpha carbon that is not occupied by atoms
     */
    
    mstreal fracUnoccupiedVolume(Residue* R) {
        if (protein == nullptr) MstUtils::error("Structure must be set","volumeCalculator::fracUnoccupiedVolume");
        return 1-occupiedVolume(R)/sphereVol(qR);
    }
    
    /**
     Same as above, just not normalized
     */
    mstreal occupiedVolume(Residue* R);
    
    mstreal querySphereVol() {return sphereVol(qR);}
    
protected:
    /**
     The good ole equation for the volume of a sphere
     */
    mstreal sphereVol(mstreal radius) {
        return (4*3.1415*radius*radius*radius)/3;
    }
    
    /**
     Finds the volume of a cap/dome of a sphere. This is the section of a sphere cut off by a plane. See  https://en.wikipedia.org/wiki/Spherical_cap
     
     The integral of (r^2-x^2)dx with limits [-r,r] gives the volume of a whole sphere, but when the limits are changed to [r-h,r] (where h is
     the distance from the surface of the sphere to the plane that intersects it) it gives the volume of the cap.
     
     @param r the radius of the whole sphere
     @param cH cap height: the distance from the surface of the sphere to the plane that intersects it
     @return the volume of the cap/dome
     */
    mstreal sphereCapVol(mstreal r, mstreal cH) {
        return (3.1415*cH*cH*(2*r-cH))/3;
    }
    
    /**
     Finds the point at which two spheres intersect.
          
     Two spheres are defined as follows (with offsets so they can be centered anywhere):

        (x-x_1)^2 + (y-y_1)^2 + (z-z_1)^2 = r_1^2 and (x-x_2)^2 + (y-y_2)^2 + (z-z_2)^2 = r_2^2
     
     First we set sphere 1 on the origin and sphere 2 with y and z at zero.
     
        x^2 + y^2 + z^2 = r_1^2 and (x-x_2)^2 + y^2 + z^2 = r_2^2

     We solve these for y and set them equal
     
        r_1^2 - x^2 = r_2^2 - (x-x_2)^2
     
     Now we solve for x, the distance at which they intersect
     
        x = (r_1^2 + x_2^2 - r_2^2)/2*x_2
     
     @param r1 the radius of the sphere on the origin
     @param r2 the radius of the sphere off the origin
     @param d the distance between the spheres
     @return the distance at which the spheres intersect
     */
    mstreal sphereIntersectDistance(mstreal r1, mstreal r2, mstreal d){
        return (r1*r1+d*d-r2*r2)/(2*d);
    }
    
    /**
     Finds the volume of two intersecting spheres.
     
     The volume of two intersecting spheres can be found by summing the volume of the two "caps" formed on either side of the intersection
     
     see the image here https://math.stackexchange.com/a/1561061
     
     @param r1 the radius of the sphere on the origin
     @param r2 the radius of the sphere off the origin
     @param d the distance between the spheres
     @return the volume of the intersection of the spheres
     */
    mstreal intersectingSphereVol(mstreal r1, mstreal r2, mstreal d) {
        mstreal intersectDistance = sphereIntersectDistance(r1,r2,d);
//        cout << "radii: " << r1 << " " << r2 << endl;
//        cout << "distance: " << d << endl;
//        cout << "intersectDistance: " << intersectDistance << endl;
//        cout << "cap height: " << r1-intersectDistance << " " << intersectDistance-(d-r2) << endl;
//        cout << "sphereCapVol: " << sphereCapVol(r1,r1-intersectDistance) << " " << sphereCapVol(r2,intersectDistance-(d-r2)) << endl;
        return sphereCapVol(r1,r1-intersectDistance) + sphereCapVol(r2,intersectDistance-(d-r2));
    }
    
private:
    Structure* protein = nullptr;
    vector<Atom*> atoms;
    ProximitySearch ps;
    
    mstreal qR = 12.0; // the radius of the sphere we will search around the alpha carbon
};

// Implements the Needleman-Wunsch algorithm for finding the optimal global alignment of backbone atoms
class peptideAlignmentNW {
public:
    peptideAlignmentNW(mstreal _gap_penalty = -20) : gap_penalty(_gap_penalty) {};
    
    // This method assumes both structures have a single chain representing the peptide and only
    // backbone atoms.
    mstreal findOptimalAlignment(Structure* _S1, Structure* _S2);

    static mstreal getDistance(vector<Residue*> R1, vector<Residue*> R2);

    static mstreal similarity(Residue* R1, Residue* R2);
private:
    mstreal gap_penalty;
    Structure* S1;
    Structure* S2;
    vector<Residue*> R1;
    vector<Residue*> R2;
    vector<vector<mstreal>> score_table;
    vector<vector<pair<int,int>>> traceback_table;
    list<pair<int,int>> alignment;
    vector<Residue*> R1_aligned;
    vector<Residue*> R2_aligned;
    mstreal best_alignment_distance;
    bool ready = false;
    bool verbose = true;
};

// Miscellaneous useful functions

/**
 Equivalent to Python's str.split()
 */
vector<string> splitString(string s, const string &delim);

/**
 Equivalent to Python's str.join(), but the components can be any
 type that can be fed into a stringstream (e.g. using <<)
 */
template<typename Streamable>
string joinString(const vector<Streamable> &components, const string &delim) {
    stringstream ss;
    for (int i = 0; i < components.size(); ++i) {
        ss << components[i];
        if (i < components.size() - 1) {
            ss << delim;
        }
    }
    return ss.str();
}

/**
 Determines the target residue from which the given seed came. Assumes the file
 name is of the form "XXXX_XABC_match_XX_seed_X.pdb", and extracts the number
 "ABC".
 
 @param seedName the file name of the seed
 @return the target residue index for this seed, or -1 if one was not found
 */
int getTargetResidueIndex(string seedName);

/**
 * Determines the target residue code, including the chain ID and the number,
 * assuming the file name is of the form "[TARGETNAME]-[CHAINID][NUMBER]-etc.pdb".
 */
pair<string, int> getTargetResidueCode(string seedName);


// =========== Utilities inherited from structgen, added 6/15/2020

// makes an ordered pair from two elements
template <class T1, class T2>
pair<T1, T2> orderedPair(T1 a, T2 b) {
    if (!(a < b)) { // greater than
        return make_pair(b, a);
    }
    return make_pair(a, b); //less than or equals to
    
}

// hash function for a pair - so that it can be used in an unordered set or map
struct pair_hash {
    template <class T1, class T2>
    size_t operator () (const pair<T1,T2> &p) const {
        auto h1 = hash<T1> {} (p.first);
        auto h2 = hash<T2> {} (p.second);
        return h1 ^ h2;
    }
};

// tests whether atom is backbone
bool isBackboneHeavy(Atom& atom);

bool hasBackbone(Residue& res, bool requireOxygen = true);

double bondDistance(Residue& res1, Residue& res2);

// is res1 close enough to be bonded to res2?
// order matters here
bool isBonded(Residue& res1, Residue& res2, double maxPepBond = 2.0);

// removes the side chains from a structure
void removeSideChains(Structure& structure);

// returns a new AtomPointerVector without the side chains
AtomPointerVector backboneAtoms(AtomPointerVector& apv);

// removes the hydrogens from a structure
void removeHydrogens(Structure& s);

// an alphabet for naming chains
string generateChainID(int i);

// combines a vector of structures into a single structure where each former structure is now a chain
Structure combineStructures(const vector<Structure>& structs, bool split = false, bool renumber = true);

// splits a structure into a vector of structures based on the chains
vector<Structure> splitByChains(const Structure& s);

string residueID(Residue& res, string sep = "", bool full = true);

// writes residues of a vector out
// perhaps make a write for a vector of pointers (template function)
void writeResidues(ostream& os, const vector<Residue*>& residues);

// writes residues from a vector out
void writeResidues(ostream& os, const vector<vector<Residue*> >& residues);

vector<Residue*> sortByResidueIndex(vector<Residue*>& residues);

vector<Residue*> sortByResidueIndex(set<Residue*>& residues);

map<Residue*, int> residueIndexMapping(vector<Residue*>& residues, bool fromResidue = true);

AtomPointerVector residuesToAtoms(vector<Residue*>& residues);

Matrix seqFreqs(const vector<Sequence>& seqs, double psuedocount = 0.0, bool normalize = true, bool nonStd = false);

Matrix seqFreqs(const vector<Sequence>& seqs, Matrix& pseudocounts, bool normalize = true);

// TODO
// remove chains below a certain size (e.g. 2 residues)
// option for dealing non standard AAs
// incorporate TargetStruct redundancy here
void cleanStructure(Structure& s, Structure& cleaned, bool reassignChains = true, bool renumber = true, int minSegLen = 1);

Structure cleanStructure(Structure& s, bool reassignChains = true, bool renumber = true, int minSegLen = 1);

// vector of chain lengths for a structure
vector<int> chainLengths(Structure& s);

// calulate density
double calculateDensity(Structure& s);

// calculate pairwise RMSD
void pairwiseRMSD(vector<AtomPointerVector>& apvs, unordered_map<int, unordered_map<int, double> >& dists, double rmsdCut);

set<pair<Atom*, Atom*> > findClashes(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& qAPV, double ratio = 0.75);

// simple clash function
// does not take into account atom sizes, ect
// no ability to exclude certain atoms or residues
bool isClash(ProximitySearch& ps, AtomPointerVector& queryAPV, double clashDist = 2.0, int maxNumClashes = 0);

// determines whether there a clash betwen an AtomPointerVector (apv) and another structure represented by a ProximitySearch object (ps)
// ps is a ProximitySearchObject
// apv is an AtomPointerVector
// psExclude is the set of atom indices in psAPV (and ps) to exlude from clashes
// ratio of sum of vdw radii (default = 0.7)
// maxNumClasshes is the maximum number of allowed clashes between queryAPV and psAPV (default = 0)
// /home/ironfs/scratch/grigoryanlab/cmack2357/minCover/lists
bool isClash(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& queryAPV, set<pair<Residue*, Residue*> >& exclude, double ratio = 0.75, int maxNumClashes = 0);

// calls the above function without excluded residues
bool isClash(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& queryAPV, double ratio = 0.75, int maxNumClashes = 0);

// Added by sebastian 9/29/20
// Moved the logic of isClash() into this function, which can be called directly if the user wants
// to know how many clashes are counted.
// if maxNumClashes is < 0, all clashes will be counted.
// Note: This function has been modified such that it reports the number of unique pairs of
// residues with a clash
int numClash(ProximitySearch& ps, AtomPointerVector& psAPV, AtomPointerVector& queryAPV, set<pair<Residue*, Residue*> >& exclude, double ratio = 0.75, int maxNumClashes = 0);

/*
 After looking through the function implementation, it's not clear to me that this function is
 necessary. It appears that the same functionality can be achieved by calling with the above
 function with the ProximitySearch, psAPV, and queryAPV all defined using the same atoms.
 */
//// Added by venkats 1/28/20
//// Same as above, but operates on a single structure
//bool isClashSingleStructure(ProximitySearch& ps, AtomPointerVector& psAPV, double ratio = 0.7, int maxNumClashes = 0);
//
//// Added by sebastian 9/29/20
//// Same as numClash, but for single structure
//// if maxNumClashes is < 0, all clashes will be counted
//int numClashSingleStructure(ProximitySearch& ps, AtomPointerVector& psAPV, double ratio = 0.7, int maxNumClashes = 0);

vector<int> effectiveSegLengths(Structure& s, vector<Residue*> subset, int minGapCont, int minSegLen = 3);

Structure residuesToStructure(vector<Residue*>& residues, double maxPeptideBond = 2.0, int startChainIdx = 0);

// Contact list utils

// writes the contents of a contactList
void writeContactList(ostream& os, contactList& cl);

// get the backbone backbone contacts for a structure
// queryResidues should be a subset of the residues in s
contactList backboneContacts(Structure& s, double dist, vector<Residue*> queryResidues = {});


contactList vdwBackboneContacts(Structure& s, double lb, double ub, vector<Residue*> queryResidues = {});

// union of contact list
// does not allow duplicates
// orders the contacts - perhaps give an option to not do this
contactList contactListUnion(const vector<contactList>& contLists, bool order = true);

#endif /* utilities_h */
