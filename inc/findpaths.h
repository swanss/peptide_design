//
//  findpaths.h
//  dtermen
//
//  Created by Venkatesh Sivaraman on 2/5/19.
//

#ifndef findpaths_h
#define findpaths_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <list>

#include "msttypes.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "fusehelper.h"
#include "seedgraph.h"
#include "seedmap.h"
#include "ranklist.h"

class PathFinder {
public:
    PathFinder(SeedGraph *g, SeedGraphMap<mstreal> *scores): graph(g), scores(scores) {
        progressPrintInterval = (int) g->residueSize() / 100;
    };
    
    vector<pair<vector<Residue *>, mstreal>> findPaths(int numPaths);
    
private:
    SeedGraph *graph;
    SeedGraphMap<mstreal> *scores;
    int progressPrintInterval;
    unordered_map<Residue *, RankList<Residue *>> pathPointers;
    unordered_set<Residue *> temporaryMarks;
    
    /**
     Computes the k best-weight paths that contain the given residue.
     */
    void computeBestPaths(Residue *residue, int numPaths);
    
    /**
     Constructs a path corresponding to the given rank index
     in the given rank list.
     
     For example, given the following graph:
     
     @code
     A --> C --> D
         /'  \,
        B     E
     @endcode
     
     and rank lists (null corresponds to path termination):
     
     @code
     A: C, C, null
     B: C, C, null
     C: D, E, null
     D: null
     E: null
     Overall: A, B, A
     @endcode
     
     Calling constructPath on the overall rank list with a
     rank index of 0 would yield the path A, C, D; with 1,
     B, C, D; and with 2, A, C, E.
     
     @param rankList a rank list of residues
     @param rankIndex the index in the rank list to construct
            the path for
     @return a path representing the rankIndex'th best-scoring
             path according to this rank list
     */
    vector<Residue *> constructPath(RankList<Residue *> &rankList, int rankIndex);
};

#endif /* findpaths_h */
