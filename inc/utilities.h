//
//  utilities.h
//  TPD_target
//
//  Created by Sebastian Swanson on 6/14/19.
//

#ifndef utilities_h
#define utilities_h

#include "msttypes.h"
#include "mstcondeg.h"
#include <unordered_map>

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

/* Borrowed from Venkat */


class StructuresBinaryFile {
public:
  StructuresBinaryFile(string filePath, bool read = true): _filePath(filePath), readMode(read) {
    openFileStream(filePath);
  };
  
  StructuresBinaryFile(const StructuresBinaryFile &other): readMode(true), _filePath(other._filePath), _filePositions(other._filePositions), _structureNames(other._structureNames) {
    MstUtils::assert(other.readMode, "Copying write-only binary file not supported");
    cout << "Opening file stream for copy, with " << _structureNames.size() << " loaded structure names" << endl;
    openFileStream(_filePath);
  }
  
  ~StructuresBinaryFile() {
    fs.close();
  }
  
  //Get structures
  Structure* next();
  Structure* getStructureNamed(string name);
  
  //Get properties
  size_t structureCount() {
    if (_filePositions.empty())
      scanFilePositions();
    return _filePositions.size();
  }

  //Navigate the file
  bool hasNext();
  void skip();
  void jumpToStructureIndex(int idx);
  void reset();
  
  //Setters
  void appendStructure(Structure *s);
  //not sure about this
  //  void insertStructureNames(vector<string> &names) {
  //    if (_filePositions.empty())
  //      scanFilePositions();
  //    names.insert(names.end(), _structureNames.begin(), _structureNames.end());
  //  }
  
protected:
  void scanFilePositions();
  void openFileStream(string filePath);

private:
  string _filePath;
  bool readMode;
  fstream fs;
  unordered_map<string, long> _filePositions; // Positions of each structure in the file
  vector<string> _structureNames;
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
};

#endif /* utilities_h */
