#include <stdio.h>
#include "findfuseable.h"
#include "mstsystem_exts.h"
#include "structure_iter.h"
#include "seedscore.h"
#include "secondarystructure.h"
#include "clustertree.h"
#include "clusterutils.h"
#include "utilities.h"
#include <unordered_set>
#include "msttypes.h"

using namespace std;

static vector<string> bba = { "CA", "O", "N", "C" };


/**
 Writes a line of the point/line cloud file.
 
 @param s the structure whose coordinates to write
 @param out the file stream to write to
 @param chainID the ID of the chain to write, or "" if writing all chains
 @param atomMode the mode specifying which atoms to write ('ca' or 'bb')
 @param isPointCloud if true, writes point cloud format, otherwise writes line cloud format
 @param colorFn a callback function that writes color information. Takes three arguments: the
    output file stream, a secondary structure classification, and the current atom.
 */
template <typename ColorFunction>
void const writeStructureCoordinates(const Structure &s, secondaryStructureClassifier &ssClassifier, ofstream &out, string chainID, string atomMode, bool isPointCloud, ColorFunction colorFn) {
    for (int i = 0; i < s.chainSize(); i++) {
        Chain &chain = s.getChain(i);
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
                colorFn(out, ssID, a);
                if (isPointCloud || res->getResidueIndexInChain() == chain.residueSize() - 1)
                    out << endl;
                else
                    out << ",";
            }
        }
    }
}

Structure getSeedSegment(Structure* S, string seed_chain = "0") {
    Chain* C = S->getChainByID(seed_chain);
    return Structure(*C);
}

int main (int argc, char *argv[]) {
    
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Samples seeds for display. Either builds a point cloud file suitable for opening with SeedViewer or exports as PDB files.");
    opts.addOption("structures", "Path to a binary file containing structures to plot", true);
    opts.addOption("clusters", "Path to a cluster tree text file generated with a single fragment fetcher", false);
    opts.addOption("segLen", "Length of the segments used in the cluster tree clustering (default = 3).", false);
    opts.addOption("chain", "Chain ID to use (default is to add all chains)", false);
    opts.addOption("atoms", "'ca' to write only alpha carbons (default); 'bb' to write all backbone atoms", false);
    opts.addOption("color", "'uniform' to write the same color for all atoms (default); 'atom' to color by atom; 'ss' to color by secondary structure. If --clusters is provided, writes out the cluster index path and ignores this parameter.", false);
    opts.addOption("out", "Path to name output file (.pcloud, .lcloud). If not provided, exports .pdb files", false);
    opts.addOption("sample", "Number of structures to subsample (default 1000)", false);
    opts.setOptions(argc, argv);

    string binaryFilePath = opts.getString("structures");
    string outputPath = opts.getString("out","");
    string colorMode = opts.getString("color", "atom");
    string atomMode = opts.getString("atoms", "ca");
    string chainID = opts.getString("chain", "");
    int sample = opts.getReal("sample", 1000);
        
    unordered_map<string, int> atomColors;
    atomColors["CA"] = 0;
    atomColors["C"] = 1;
    atomColors["O"] = 2;
    atomColors["N"] = 3;

    // Check output path file extension
    string outFileExt = outputPath.substr(outputPath.find_last_of(".") + 1);
    bool pdb = false;
    if (outFileExt != "pcloud" && outFileExt != "lcloud") pdb = true;
    if (outFileExt == "lcloud" && atomMode == "bb") MstUtils::error("lcloud and bb atom mode not supported");

//    ofstream out(outputPath, ios::out);
//    if (!out.is_open())
//        MstUtils::error("Could not open file stream");
//    out << std::fixed << std::setprecision(4) << endl;

    StructuresBinaryFile binaryFile(binaryFilePath);
    mstreal rate = min(sample/mstreal(binaryFile.structureCount()),1.0);
    secondaryStructureClassifier ssClassifier;

    if (opts.isGiven("clusters")) {
//        // We want to render the cluster tree by enumerating the cluster nodes, then listing
//        // the coordinates as well as the address for each node. For example, the second child
//        // of the first child of the root node would be listed as 1:2.
//        SingleFragmentFetcher fetcher(&binaryFile, opts.getInt("segLen", 3), chainID);
//        ClusterTree tree(&fetcher, 4, true); // shared coordinates
//        cout << "Reading cluster tree..." << endl;
//        tree.read(opts.getString("clusters"));
//
//        auto root = tree.getRoot();
//        int treeSize = root->subtreeSize();
//
//        ClusterIterator it(root);
//        int nodeIndex = 0;
//        while (it.hasNext()) {
//            auto node = it.next();
//            if (nodeIndex++ % 1000 == 0)
//               cout << "Cluster element " << nodeIndex << endl;
//
//            if (rand() / (float)RAND_MAX >= rate) {
//                continue;
//            }
//
//            Structure s = fetcher.getResultStructure(node->getItem());
//            string address = joinString(it.currentAddress(), ":");
//
//            writeStructureCoordinates(s, ssClassifier, out, "", atomMode, outFileExt == "pcloud", [&](ofstream &out, string ssID, Atom *a) {
//                out << address;
//            });
//        }
//
    } else {
        // The simple case is to enumerate all the structures in the binary file, and write out
        // their coordinates with the appropriate color mode.
        int structureIdx = 0;
        while (binaryFile.hasNext()) {
            if (structureIdx++ % 1000 == 0)
               cout << "Structure " << structureIdx << endl;

            if (rand() / (float)RAND_MAX >= rate) {
                binaryFile.skip();
                continue;
            }

            Structure *s = binaryFile.next();
            Structure S_seed = getSeedSegment(s);
            if (pdb) {
                S_seed.writePDB(s->getName()+".pdb");
            }
//            else {
//                writeStructureCoordinates(S_seed, ssClassifier, out, chainID, atomMode, outFileExt == "pcloud", [&](ofstream &out, string ssID, Atom *a) {
//                    if (colorMode == "uniform") {
//                        out << 0;
//                    } else if (colorMode == "atom") {
//                        out << atomColors[a->getName()];
//                    } else if (colorMode == "ss") {
//                        out << ssID;
//                    }
//                });
//            }
            delete s;
        }
    }
        
//    out.close();
    return 0;
}
