#include <stdio.h>
#include <fstream>
#include "mstoptions.h"
#include "msttypes.h"
#include "structure_iter.h"
#include "utilities.h"
#include "mstsystem_exts.h"

bool interactionSort(const tuple<int, int>& a, 
               const tuple<int, int>& b) {
    return (get<1>(a) > get<1>(b));
}

void randomSelection(int numStructures, int totalStructures, vector<int>& structureInds) {
    for (int i = 0; i < numStructures; i++) {
        bool newInt = false;
        while (!newInt) {
            int randInt = MstUtils::randInt(0,totalStructures-1);
            if (find(structureInds.begin(), structureInds.end(), randInt) == structureInds.end()) {
                newInt = true;
                structureInds.push_back(randInt);
            }
        }
    }
}

int main(int argc, char* argv[]) {
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Selects structures either randomly or according to their number of contacts and in a way that maximizes structural diversity.");
    opts.addOption("mode", "Indicator for random (0) or procedural (1) selection.");
    opts.addOption("num", "Number of structures to select.");
    opts.addOption("total", "Total number of structures to select from.");
    opts.addOption("list", "Path to file listing each of the path structures.");
    opts.addOption("distanceList", "Path to file listing the relative distances between each structure.");
    opts.addOption("interactionList", "Path to file listing the number of target interactions for each structure.");
    opts.addOption("name", "PDB identifier of input structure.");
    opts.addOption("threshold", "Angstrom threshold on RMSD distances between chosen structures.");
    opts.addOption("minContacts", "Minimum number of contacts for chosen structures.");
    opts.addOption("out", "Output directory.");
    opts.setOptions(argc, argv);

    int numStructures = opts.getInt("num");
    vector<int> structureInds;
    int totalStructures = opts.getInt("total");
    // random selection
    if (opts.getInt("mode") == 0) {
        randomSelection(numStructures, totalStructures, structureInds);        
    } else { // procedural selection
        // sort structures by number of interactions
        fstream fin;
        string interactionsPath = opts.getString("interactionList");
        vector< tuple<int, int> > interactions;
        vector<string> structureData = MstUtils::fileToArray(interactionsPath);
        string header = structureData[0];
        int numInteractionsPos = 0;
        while (true) {
            string curHeader = MstUtils::nextToken(header, ",", true);
            if (curHeader.compare("num_interface_contacts") == 0) {
                break;
            }
            numInteractionsPos++;
        }
        structureData.erase(structureData.begin());
        vector<string> rowData;
        int count = 0;
        string data;
        for (string line : structureData) {
            for (int i = 0; i <= numInteractionsPos; i++) {
                data = MstUtils::nextToken(line, ",", true);
            }
            interactions.push_back(std::make_tuple(count, MstUtils::toInt(data)));
            count++;
        }
        sort(interactions.begin(), interactions.end(), interactionSort);

        // read in distance matrix
        string distancesPath = opts.getString("distanceList");
        vector<string> distancesData = MstUtils::fileToArray(distancesPath);
        map<int, int> headerData;
        string headerName;
        string headerLine = distancesData[0];
        string headerInfo;
        count = 0;
        while (!headerLine.empty()) {
            headerInfo.clear();
            headerName = MstUtils::nextToken(headerLine, ",", true);
            while (!headerName.empty()) {
                headerInfo = MstUtils::nextToken(headerName, "_", true);
            }
            headerData.insert(pair<int, int>(MstUtils::toInt( MstUtils::nextToken(headerInfo, ".", true)), count));
            count++;
        }
        distancesData.erase(distancesData.begin());
        vector< vector<mstreal> > distances;
        vector<mstreal> rowDistances;
        for (string line : distancesData) {
            rowDistances.clear();
            MstUtils::nextToken(line, ",", true);
            while (!line.empty()) {
                rowDistances.push_back(MstUtils::toReal(MstUtils::nextToken(line, ",", true)));
            }
            distances.push_back(rowDistances);
        }

        // greedily choose interactions that are sufficiently different from already chosen structures
        int numRemaining = numStructures;
        int searchedNum = 0;
        int curStructureInd;
        int minContacts = opts.getInt("minContacts");
        mstreal threshold = opts.getReal("threshold");
        while (numRemaining > 0) {
            curStructureInd = get<0>(interactions[searchedNum]);
            if (get<1>(interactions[searchedNum]) < minContacts) {
                cout << "No more structures with enough contacts after " << (numStructures - numRemaining) << " have been found." << endl;
                break;
            }
            mstreal curStructureDistance = INT_MAX;
            for (int chosenStructureInd : structureInds) {
                curStructureDistance = min(curStructureDistance, distances[headerData[curStructureInd]][headerData[chosenStructureInd]]);
            }
            if ((curStructureDistance > threshold) || (curStructureDistance == INT_MAX)) {
                structureInds.push_back(curStructureInd);
                numRemaining--;
            }
            if ((searchedNum == totalStructures) && (numRemaining > 0)) {
                cout << "No more sufficiently different structures after " << (numStructures - numRemaining) << " have been found." << endl;
                //cout << "No more sufficiently different structures after " << (numStructures - numRemaining) << " have been found. Resorting to random selection for the remaining " << numRemaining << " structures." << endl;
                //randomSelection(numRemaining, totalStructures, structureInds);
                numRemaining = 0;
            }
            searchedNum++;
        }
    }

    // write names of chosen structures to structures.list
    string structureName = opts.getString("name");
    string outputFile = MstSystemExtension::join(opts.getString("out"), "structures.list");
    if (MstSys::fileExists(outputFile)) MstSys::crm(outputFile);
    string outputDir = MstSystemExtension::join(opts.getString("out"), "structures/");
    MstSys::cmkdir(outputDir);
    fstream out;
    MstUtils::openFile(out, outputFile, ios_base::out);
    for (int structureInd : structureInds) {
        string outputName = structureName + "_fused-path_" + MstUtils::toString(structureInd) + ".pdb";
        out << outputName << endl;
    }
    out.close();
        
    cout << "Done" << endl;

    return 0;
}
