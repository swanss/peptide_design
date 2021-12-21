
#ifndef MAKECLUSTERTESTSET_H
#define MAKECLUSTERTESTSET_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <utility>
#include <sys/types.h>
#include <unistd.h>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <climits>
#include <cfloat>
#include <stdio.h>


#include "msttypes.h"
#include "msttransforms.h"
#include "mstrotlib.h"
#include "mstcondeg.h"
#include "mstfasst.h"
#include "mstlinalg.h"
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "mstfuser.h"
#include "structure_iter.h"
#include "seedgraph.h"
#include "seedscore.h"
#include "utilities.h"
#include "params.h"
#include "dtermen.h"

using namespace std;
using namespace MST;


class clusterTestSet {
    public:
        clusterTestSet(string& _structuresInputPath, string& _structuresList, Structure* _target, Structure* _peptide, string& configPath, int _nStructures): config(configPath) {
            // StructureCache* structuresCache = new StructureCache(_structuresInputPath);
            // structuresCache->preloadFromPDBList(_structuresList);
            structuresPath = _structuresInputPath;
            structuresList = _structuresList;
            nStructures = _nStructures;
            target = _target;
            peptide = _peptide;
            rl = new RotamerLibrary(config.getRL());
        }

        void selectStructures(int negType);

        void writeStructures(string& outputPath);

        int getNumTruePositives() {
            return truePositives.size();
        }

        int getNumTrueNegatives() {
            return trueNegatives.size();
        }

        int getNumUnkownSamples() {
            return unkownSamples.size();
        }

        mstreal getDistLimit() {
            return distLimit;
        }

        void resetCombinedStructure();

        bool prepareCombinedStructure(Structure *seed);

        Residue* getTargetResidue(Residue* targetCopyRes) {
            /*
            since the seed chain is always added as the final chain, the index of the residues
            in the previous chains (i.e. the target) should not change
            */
            return &target->getResidue(targetCopyRes->getResidueIndex());
        }

    private:
//        StructureCache* structuresCache;
        string structuresPath;
        string structuresList;
        Structure* target;
        Structure* peptide;
        int nStructures;
        set<string> truePositives;
        set<string> trueNegatives;
        set<string> unkownSamples;
        configFile config;
        RotamerLibrary *rl;
        mstreal distLimit;
};


#endif
