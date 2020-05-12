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
    op.addOption("Nc","The total number of independent MC trajectories (default 1000)");
    op.addOption("Ni","The total number of iterations in a single MC trajectory (default 100,000)");
    op.addOption("base","The base name of the list of solutions",true);
    op.setOptions(argc,argv);
    
    //import parameters
    int Nc = op.getInt("Nc",1000);
    int Ni = op.getInt("Ni",100000);
    
    //import energy table
    EnergyTable E;
    string etabFile = op.getString("etab");
    E.readFromFile(etabFile);
    
    //open file to write solutions
    
    // set up output file
    string out_file = op.getString("base") + ".sols.txt";
    fstream solutions;
    MstUtils::openFile(solutions, out_file, fstream::out);
    
    /* mc has an option run multiple cyles, but as far as I'm aware, only provides the top scoring
     sequence as the solution. As a workaround, I'm going to just run it Nc times and store the sequence
     each time */
    Sequence bestSeq;
    for (int i = 0; i < Nc; i++) {
        vector<int> bestSol = E.mc(1, Ni, 1.0);
        mstreal lowE = E.scoreSolution(bestSol);
        cout << "lowest energy found is " << lowE << endl;
        bestSeq = E.solutionToSequence(bestSol);
        cout << "lowest-energy sequence: " << bestSeq.toString() << endl;
        cout << "mean energy is " << E.meanEnergy() << endl;
        cout << "estimated energy standard deviation is " << E.energyStdEst() << endl;
        
        //store sequence
        solutions << bestSeq.toString() << endl;
    }
};
