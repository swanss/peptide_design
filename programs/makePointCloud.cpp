//
//  makePointCloud.cpp
//  DummyTarget
//
//  Created by Venkatesh Sivaraman on 4/22/19.
//

#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "structure_iter.h"
#include "fusehelper.h"
#include "seedscore.h"
#include "secondarystructure.h"
#include <unordered_set>

using namespace std;

int main (int argc, char *argv[]) {
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Builds a point cloud file suitable for opening with SeedViewer.");
    opts.addOption("structures", "Path to a binary file containing structures to plot", true);
    opts.addOption("chain", "Chain ID to use (default is to add all chains)", false);
    opts.addOption("atoms", "'ca' to write only alpha carbons (default); 'bb' to write all backbone atoms", false);
    opts.addOption("color", "'uniform' to write the same color for all atoms (default); 'atom' to color by atom; 'ss' to color by secondary structure", false);
    opts.addOption("out", "Path to name output file (.pcloud or .lcloud)", true);
    opts.addOption("rate", "Fraction of structures to sample (default 1.0)", false);
    opts.setOptions(argc, argv);

    string binaryFilePath = opts.getString("structures");
    string outputPath = opts.getString("out");
    string colorMode = opts.getString("color", "atom");
    string atomMode = opts.getString("atoms", "ca");
    string chainID = opts.getString("chain", "");
    float rate = opts.getReal("rate", 1.0);

    static vector<string> bba = { "CA", "O", "N", "C" };
    unordered_map<string, int> atomColors;
    atomColors["CA"] = 0;
    atomColors["C"] = 1;
    atomColors["O"] = 2;
    atomColors["N"] = 3;

    // Check output path file extension
    string outFileExt = outputPath.substr(outputPath.find_last_of(".") + 1);
    if (outFileExt != "pcloud" && outFileExt != "lcloud")
        MstUtils::error("Must output a .pcloud or .lcloud file, not " + outFileExt);
    if (outFileExt == "lcloud" && atomMode == "bb")
        MstUtils::error("lcloud and bb atom mode not supported");

    ofstream out(outputPath, ios::out);
    if (!out.is_open())
        MstUtils::error("Could not open file stream");
    out << std::fixed << std::setprecision(4) << endl;

    StructuresBinaryFile binaryFile(binaryFilePath);

    secondaryStructureClassifier ssClassifier(7, 8);

    int structureIdx = 0;
    while (binaryFile.hasNext()) {
        if (structureIdx++ % 1000 == 0)
           cout << "Structure " << structureIdx << endl; 

        if (rand() / (float)RAND_MAX >= rate) {
            binaryFile.skip();
            continue;
        } 

        Structure *s = binaryFile.next();
        for (int i = 0; i < s->chainSize(); i++) {
            Chain &chain = s->getChain(i);
            if (!chainID.empty() && chain.getID() != chainID)
                continue;
            for (Residue *res: chain.getResidues()) {
                vector<Atom *> atoms;
                if (atomMode == "ca") {
                    Atom *ca = res->findAtom("CA");
                    if (ca)
                        atoms.push_back(ca);
                } else if (atomMode == "bb") {
                    for (string code: bba) {
                        Atom *a = res->findAtom(code);
                        if (a)
                            atoms.push_back(a);
                    }
                } else {
                    MstUtils::error("Unrecognized atom mode " + atomMode);
                }

                string ssID = ssClassifier.classification2ColorID(ssClassifier.classifyResidue(res));
                for (Atom *a: atoms) {
                    // Add it to the file
                    out << a->getX() << "," << a->getY() << "," << a->getZ() << ",";
                    if (colorMode == "uniform") {
                        out << 0;
                    } else if (colorMode == "atom") {
                        out << atomColors[a->getName()];                        
                    } else if (colorMode == "ss") {
                        out << ssID;
                    }
                    if (outFileExt == "pcloud" || res->getResidueIndexInChain() == chain.residueSize() - 1)
                        out << endl;
                    else
                        out << ",";
                }
            }
        }

        delete s;
    }
        
    out.close();
    return 0;
}
