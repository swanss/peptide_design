#ifndef fileformats_h
#define fileformats_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <set>
#include <map>
#include <list>
#include <unordered_map>
#include <fstream>

#include "msttypes.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "structure_iter.h"

using namespace std;
using namespace MST;

/**
 Type that encapsulates a CSV file indicating paths for seeds. At minimum, the
 file should contain the path for each seed and a string indicating which chains
 can be used as seeds (e.g. "AB" for both chains A and B). It may also contain
 one or more scores for these seeds.
 
 The CSV file should be formatted with a header line, where the first column is
 called "path" and the second column is called "chids".
 */
class SeedListFile {
public:
    SeedListFile(string path): _path(path) {}
    ~SeedListFile();
    
    /**
     Reads the file.
     
     @param pdbPrefix a path to prepend to each path in the file for return
     @return a list of seed paths from the file, and a list of the chain IDs
     corresponding to each seed path
     */
    pair<vector<string>, vector<string>> read(string pdbPrefix = "");
    /**
     @return the names of each metadata field in the file
     */
    vector<string> getMetadataFields() { return _metadataNames; }
    /**
     @return the stored metadata for the given path
     */
    vector<string> getMetadata(string path) { return metadata[path]; }
    /**
     Sets the metadata field names for writing in the CSV file.
     
     @param metadataNames the new field names
     */
    void setMetadataFields(vector<string> metadataNames) { _metadataNames = metadataNames; }
    /**
     Writes the given seed path to the file, with the optional scores provided.
     
     @param seed the seed path to write
     @param chainIDs a concatenated string of the seed chain IDs
     @param metadata a list of metadata items to provide in the file (must have
     previously provided a set of headers of the same length)
     */
    void write(string seed, string chainIDs, vector<string> metadata = vector<string>());
private:
    vector<string> _metadataNames;
    vector<string> listedPaths;
    vector<string> listedChainIDs;
    map<string, vector<string>> metadata; // keys paths to metadata fields
    string _path;
    ifstream *readstream = nullptr;
    ofstream *writestream = nullptr;
};

struct FuseCandidate {
    string file1;
    string chain1;
    string file2;
    string chain2;
    Structure *structure1 = nullptr;
    Structure *structure2 = nullptr;
    
    int overlapPosition1;
    int overlapPosition2;
    int overlapSize;
    float rmsd;
    
    /**
     Sets the properties of the FuseCandidate based on the given structure.
     @param structure the structure to use as the first structure in the candidate
     @param chainID the ID of the seed chain within the structure
     */
    void setStructure1(Structure *structure, string chainID) {
        structure1 = structure;
        file1 = structure->getName();
        chain1 = chainID;
    }
    /**
     Sets the properties of the FuseCandidate based on the given structure.
     @param structure the structure to use as the second structure in the candidate
     @param chainID the ID of the seed chain within the structure
     */
    void setStructure2(Structure *structure, string chainID) {
        structure2 = structure;
        file2 = structure->getName();
        chain2 = chainID;
    }
    
    /**
     Loads the structure1 and structure2 properties of the fuse candidate from
     the candidate's metadata. The user is responsible for freeing the structures
     using freeStructures().
     @param pathPrefix a path to prepend to the file1 and file2 attributes to
     find the PDB path
     */
    void loadStructures(string pathPrefix = "") {
        structure1 = new Structure(pathPrefix + file1);
        structure2 = new Structure(pathPrefix + file2);
    }
    
    /**
     Deletes the structures associated with this fuse candidate.
     */
    void freeStructures() {
        delete structure1;
        delete structure2;
    }
};

/**
 Type that encapsulates a CSV file listing fuse candidates, allowing for read
 and write operations. Note: A FuseCandidateFile can be used for reading or
 writing, but not both.
 */
class FuseCandidateFile {
public:
    FuseCandidateFile(string path, bool shouldAppend = false): _path(path), _shouldAppend(shouldAppend) {}
    ~FuseCandidateFile();
    
    /**
     Reads up to numLines lines from the file. May return fewer if the end of
     the file is reached.
     
     @param numLines the number of lines to read
     @return a list of fuse candidates from the file, with up to numLines entries
     */
    vector<FuseCandidate> read(int numLines);
    /**
     Writes the given list of fuse candidates to the file.
     
     @param candidates the list of candidates to write
     @param parentPath the path to write other paths relative to
     */
    void write(vector<FuseCandidate> candidates, string parentPath = "");
    /**
     Writes the given fuse candidate to the file.
     
     @param candidate the candidate to write
     @param parentPath the path to write other paths relative to
     */
    void write(FuseCandidate candidate, string parentPath = "");
private:
    bool _isFileDone;
    string _path;
    bool _shouldAppend;
    ifstream *readstream = nullptr;
    ofstream *writestream = nullptr;
};

#endif /* fileformats_h */
