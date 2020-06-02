//
//  fuser.h
//  dtermen
//
//  Created by Venkatesh Sivaraman on 1/15/19.
//

#ifndef fusehelper_h
#define fusehelper_h

#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <set>
#include <map>
#include <fstream>

#include "msttypes.h"
#include "mstfuser.h"
#include "mstcondeg.h"
#include "mstsystem_exts.h"
#include "seedutils.h"

using namespace std;
using namespace MST;

class FuseHelper {
public:
    /**
     Fuses the candidates listed in the given file, and writes the output PDB
     structures to the given directory.
     
     @param file the input file containing fuse candidates
     @param outDir the directory in which to write the outputs, with trailing '/'
     @param pdbPrefix the path string with which to prepend filenames from the
            candidate file in order to read them
     */
    void writeFuseResultsFromFile(FuseCandidateFile file, string outDir, string pdbPrefix = "");
    
    /**
     Fuses the given fuse candidate. Expects the candidate's structures to be
     already loaded.
     
     @param candidate the fuse candidate to fuse
     @return the fused structure, the chains in the new structure that
             correspond to the fused seed, and the fusion output object
     */
    tuple<Structure, string, fusionOutput> performFuse(FuseCandidate candidate);
};

#endif /* fusehelper_h */
