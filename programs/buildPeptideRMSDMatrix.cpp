#include <stdio.h>
#include "mstoptions.h"
#include "msttypes.h"
#include "structure_iter.h"
#include "utilities.h"
#include "mstsystem_exts.h"

int main(int argc, char* argv[]) {
    // Get command-line arguments
    MstOptions opts;
    opts.setTitle("Computes RMSD between peptide backbones in a parallelizable manner and builds a complete distance matrix for downstream analysis. Must be run three times to generate the final matrix 1) makePeptideBin, 2) computeRMSD, 3) buildMatrix");
    opts.addOption("list", "Path to file listing each of the path structures. If only this is provided, will run in makePeptideBin mode, write these all to single binary file");
    opts.addOption("bin", "Path to a binary file containing all of the paths. Provide for computeRMSD/buildMatrix");
    opts.addOption("numWorkers", "The number of workers. Provide for computeRMSD");
    opts.addOption("worker", "The index of this worker (from 1 to numWorkers). Provide for computeRMSD");
    opts.addOption("distanceList", "Path to file listing all of the distance#.csv files. If provided with --bin, will run in buildMatrix mode and combine all into distance matrix");
    opts.setOptions(argc, argv);

    if (opts.isGiven("list")) { // makePeptideBin mode
        cout << "Running in makePeptideBin mode..." << endl;
        string paths_bin_path = opts.getString("bin","./path_structures.bin");
        
        // Load PDBs from a file and write to StructuresBinaryFile
        StructuresBinaryFile* path_bin = new StructuresBinaryFile(paths_bin_path,false,false);
        
        string list_path = opts.getString("list");
        vector<string> structure_paths = MstUtils::fileToArray(list_path);
        for (string path : structure_paths) {
            Structure path_structure(path);
            cout << path_structure.getName() << endl;
            cout << MstSystemExtension::fileName(path_structure.getName()) << endl;
            path_structure.setName(MstSystemExtension::fileName(path_structure.getName()));
            path_bin->appendStructure(&path_structure);
        }
        delete path_bin;
    } else if (!opts.isGiven("distanceList") && opts.isGiven("bin")) { // computeRMSD mode
        cout << "Running in computeRMSD mode..." << endl;
        // Select batches of structures and compute distances
        int numWorkers = opts.getInt("numWorkers",1);
        int workerIndex = opts.getInt("worker",1);
        int batchSize = 10;
        string chainID = "0";
        string paths_bin_path = opts.getString("bin");
        BatchPairStructureIterator path_batch_iterator(paths_bin_path, workerIndex - 1, numWorkers, batchSize, chainID);
        
        map<pair<int,int>,mstreal> paths_distances;
        
        int count = 0;
        while (path_batch_iterator.hasNext()) {
            pair<vector<Structure *>, vector<Structure *>> path_batch = path_batch_iterator.next();
            
            // The indices refer to the specific batch, but we want to keep track of the index within the future distance matrix
            int i = path_batch_iterator.getFirstIndex() * batchSize;
            int j = path_batch_iterator.getSecondIndex() * batchSize;
            cout << "i,j: " << i << "," << j << endl;
            
            cout << "i_batch size: " << path_batch.first.size() << " j_batch size: " << path_batch.second.size() << endl;
            // i_batch/j_batch refer to the index within the batch.
            for (int i_batch = 0; i_batch < path_batch.first.size(); i_batch++) {
                Chain* path_1 = path_batch.first[i_batch]->getChainByID("0");
                for (int j_batch = 0; j_batch < path_batch.second.size(); j_batch++) {
                    count++;
                    //avoid unecessary comparisons (as distance is symmetric)
                    if (j_batch + j > i + i_batch) continue;
                    
                    Chain* path_2 = path_batch.second[j_batch]->getChainByID("0");
                    
                    //NW and RMSD
                    mstreal rmsd = generalUtilities::bestRMSD(path_1,path_2);
                                        
                    paths_distances[make_pair(i + i_batch,j + j_batch)] = rmsd;
                }
            }
        }
        cout << "count: " << count << endl;
        cout << "map size: " << paths_distances.size() << endl;
        
        // Write the info file containing all of the calculated distances
        fstream out;
        string info_file_path = "distances" + MstUtils::toString(workerIndex) + ".csv";
        MstUtils::openFile(out, info_file_path, fstream::out);
        
        for (auto it = paths_distances.begin(); it != paths_distances.end(); it++) {
            out << (*it).first.first << "," << (*it).first.second << "," << (*it).second << endl;
        }
        out.close();
    } else if (opts.isGiven("distanceList") && opts.isGiven("bin")) { // buildMatrix mode
        cout << "Running in buildMatrix mode..." << endl;
        // Read the info files and combine into a single distance matrix
        string paths_bin_path = opts.getString("bin");
        StructuresBinaryFile path_bin(paths_bin_path);
        long structure_count = path_bin.structureCount();
        vector<string> structure_names = path_bin.getStructureNames();
        string distance_list_path = opts.getString("distanceList");
        vector<string> distances_paths = MstUtils::fileToArray(distance_list_path);
        
        // Read through each distances#.csv file and add to map
        map<pair<int,int>,mstreal> paths_distances;
        for (string& path : distances_paths) {
            vector<string> lines = MstUtils::fileToArray(path);
            for (string& line : lines) {
                vector<string> split = MstUtils::split(line,",");
                if (split.size() != 3) MstUtils::error("Wrong number of elements in line from file: "+ path,"buildPathDistanceMatrix::main");
                int i = MstUtils::toInt(split[0]);
                int j = MstUtils::toInt(split[1]);
                mstreal distance = MstUtils::toReal(split[2]);
                paths_distances[pair<int,int>(i,j)] = distance;
            }
        }
        
        // Construct a proper distance matrix from the data
        fstream out;
        string info_file_path = "complete_distance_matrix.csv";
        MstUtils::openFile(out, info_file_path, fstream::out);
        //the first line is the header, and the first position is blank
        for (string structure_name : structure_names) {
            out << "," << structure_name;
        }
        out << endl;
        
        //now write the distances
        for (int i = 0; i < structure_count; i++) {
            out << structure_names[i] << ",";
            for (int j = 0; j < structure_count; j++) {
                mstreal distance;
                if (i >= j) {
                    distance = paths_distances[pair<int,int>(i,j)];
                }
                else {
                    // since each distance is only recorded once, find the value mirrored across the diagonal
                    distance = paths_distances[pair<int,int>(j,i)];
                }
                if (j < structure_count - 1) {
                    out << distance << ",";
                } else {
                    out << distance;
                }
            }
            out << endl;
        }
        out.close();
    } else {
        MstUtils::error("Arguments are not compatible with any of the three modes");
    }
    
    cout << "Done" << endl;
    
    return 0;
}
