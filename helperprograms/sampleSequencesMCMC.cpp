//
//  testMCSampleEtab.cpp
//  TPD_target
//
//  Created by Sebastian Swanson on 8/14/19.
//

#include "msttypes.h"
#include "mstoptions.h"
#include "dtermen.h"

using namespace std;
using namespace MST;

int main (int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("Given an etab, generates N top scoring sequences using a Monte Carlo Method and outputs the sequences as a flat file.");
    op.addOption("etab","The energy table that contains the self and pair energies",true);
    op.addOption("N","The desired number of final sequences default (1000)");
    op.addOption("Nc","The total number of independent MC trajectories (default 1)");
    op.addOption("Ni","The total number of iterations in a single MC trajectory (default 100,000)");
    op.addOption("kTi","The temperature when a cycle begins (default 1.0)");
    op.addOption("kTf","The temperature when the cycle ends. Ignored unless > 0 (default -1.0)");
    op.addOption("tempRange","If specified, will try different kTf and repeat the sampling (default false)");
    op.addOption("base","The base name of the file containing solutions",true);
    op.setOptions(argc,argv);
    
    //import parameters
    int N = op.getInt("N",1000);
    int Nc = op.getInt("Nc",1);
    int Ni = op.getInt("Ni",100000);
    mstreal kTi = op.getReal("kTi",1.0);
    mstreal kTf = op.getReal("kTf",-1.0);
    if (kTi < kTf) MstUtils::error("kTi must be >= kTf");
    bool tempRange = op.isGiven("tempRange");
    mstreal stepSize = 0.2;
    
    //import energy table
    EnergyTable E;
    string etabFile = op.getString("etab");
    E.readFromFile(etabFile);
    
    //open file to write solutions
    
    // set up output file
    string out_file = op.getString("base") + ".sols.txt";
    fstream solutions;
    MstUtils::openFile(solutions, out_file, fstream::out);
    
    // create another file with energies
    string csv_path = op.getString("base") + ".sols.csv";
    fstream csv;
    MstUtils::openFile(csv, csv_path, fstream::out);
    csv << "sequence,energy,Nc,Ni,kTi,kTf" << endl;
    
    /* mc has an option run multiple cyles, but as far as I'm aware, only provides the top scoring
     sequence as the solution. As a workaround, I'm going to just run it Nc times and store the sequence
     each time */
    Sequence bestSeq;
    vector<mstreal> kTf_vec;
    if (tempRange) kTf_vec = MstUtils::range(kTf,kTi,stepSize);
    else kTf_vec = {kTf};
    cout << "will try " << kTf_vec.size() << " different kTf" << endl;
    for (mstreal kTf_i : kTf_vec) {
        cout << "try kTf: " << kTf_i << endl;
        for (int i = 0; i < N; i++) {
            cout << "Sequence " << i << endl;
            vector<int> bestSol = E.mc(Nc, Ni, kTi, kTf_i);
            mstreal lowE = E.scoreSolution(bestSol);
            cout << "lowest energy found is " << lowE << endl;
            bestSeq = E.solutionToSequence(bestSol);
            cout << "lowest-energy sequence: " << bestSeq.toString() << endl;
            cout << "mean energy is " << E.meanEnergy() << endl;
            cout << "estimated energy standard deviation is " << E.energyStdEst() << endl;
            
            //store sequence
            solutions << bestSeq.toString() << endl;
            
            //store sequence/energy
            csv << bestSeq.toString() << "," << lowE << "," << Nc << ",";
            csv << Ni << "," << kTi << "," << kTf_i << ",";
            csv << endl;
        }
    }
};
