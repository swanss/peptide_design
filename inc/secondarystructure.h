//
//  secondarystructure.h
//  TPD_target
//
//  Created by Sebastian Swanson on 2/16/20.
//

#ifndef secondarystructure_h
#define secondarystructure_h

#include "msttypes.h"
#include "utilities.h"
#include "structure_iter.h"

using namespace MST;

class secondaryStructureClassifier;

/* --------- dihedralProbability --------- */

class dihedralProbabilities {
public:
    friend secondaryStructureClassifier;
    /*
     This class borrows the probabilities used in STRIDE to assign secondary structure
     PROTEINS: Structure, Function, and Genetics 23566-579 (1995)
     
     These values reflect the probability that a residue would be classified as part of a helix/sheet,
     given its dihedral angles. The dihedral angles are divided into 20x20 spaces (so the there are
     18x18 bins total).
     */
    enum SecStructType {HELIX,SHEET};
    
    dihedralProbabilities(SecStructType type);
    
    mstreal getResidueProb(Residue* R);
    
protected:
    int angle2Bin(mstreal angle);
    
private:
    vector<vector<mstreal>> probabilities;
};

/* --------- secondaryStructureClassifier --------- */

class secondaryStructureClassifier {
public:
    secondaryStructureClassifier(int _flanking_res = 2, int _agreement = 4) : flanking_res(_flanking_res), agreement(_agreement), helix(dihedralProbabilities::HELIX), sheet(dihedralProbabilities::SHEET) {
        stride = false;
        threshold = float(agreement)/float((2*flanking_res)+1);
    };
    
    /*
     Classifications are made by residue
     
     H - Helix
     E - Sheet
     O - Other
     
     In order to be classified as helix/sheet a residue must:
     
     1) have a non-zero probability of falling into the class (based on dihedral
     angles). Note: if it has a non-zero probability of falling into two classes,
     take the higher probability class
     
     2) at least some number (agreement) flanking res must have a non-zero probability of
     the same class.
     */
    
    void set_stride(bool _stride) {stride = _stride;};
    
    string classifyResidue(Residue* R);
    string classifyChain(Chain* C);
    string classifyStruct(const Structure& S);
    string classifyResInStruct(Structure* S, vector<int> res_idx);
    void assignSSToSequence(Structure* seed);
    tuple<int,int,int> getSecStructAll(const Structure& S);
    tuple<mstreal,mstreal,mstreal> getSecStructFractions(const Structure& S);
    void writeResClassificationtoFile(Chain* C, fstream& out);
    void writeResClassificationtoTSV(Structure& S, fstream& out);
    void writeCentroidtoPointFile(string bin_path, int num_final_seeds, fstream& out);
    void writeCaInfotoPointFile(Chain* C, fstream& out);
    void writeCaInfotoPointFile(string bin_path, fstream& out);
    void writeCaInfotoLineFile(Chain* C, fstream& out);
    void writeCaInfotoLineFile(string bin_path, int num_final_seeds, fstream& out);
    void writeCaInfotoLineFile(string bin_path, vector<string> seed_names, fstream& out);
    
    string classification2ColorID(string classification) {
        if (classification == "H") return "0";
        else if (classification == "E") return "1";
        else return "2";
    };
    
protected:
    
private:
    int flanking_res, agreement;
    mstreal threshold;
    bool stride; //ignore NCAA if true
    dihedralProbabilities helix, sheet;
};

#endif /* secondarystructure_h */
