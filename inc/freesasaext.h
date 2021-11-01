#ifndef freesasaext_h
#define freesasaext_h

#include "msttypes.h"
#include "mstrotlib.h"

#include "freesasa.h"

using namespace MST;

class sasaCalculator {
public:
    sasaCalculator(Structure* _parentStructure) {
        sasaDefined = false;
        params.n_threads = 1;
        setAllowedAA();
        
        setStructure(_parentStructure);
    };
    
    sasaCalculator() {
        sasaDefined = false;
        params.n_threads = 1;
        setAllowedAA();
    };

    ~sasaCalculator() {
        freesasa_structure_free(structure);
        freesasa_node_free(root_node);
    };

    void setStructure(Structure* _parentStructure);

    void computeSASA();

//    mstreal getStructureSASA();
    map<Residue*,mstreal> getResidueSASA(bool relative = false);

    map<Atom*,mstreal> getAtomSASA();

protected:
    void setAllowedAA();
    
    void setSASAStructure();

    bool sasaDefined = false; //true if sasa has been calculated for current structure
    
    void leftPadTo(string& str, const size_t num, const char paddingChar = ' ');
private:
    //mst variables
    Structure* parentStructure = nullptr;
    map<pair<string,int>,Residue*> resMap; //key is chainID and res number
    
    set<string> allowedAA;

    //freesasa variables
    freesasa_parameters params = freesasa_default_parameters;
    freesasa_structure *structure = NULL;
    const freesasa_classifier *classifier = &freesasa_naccess_classifier;
    freesasa_node *root_node = NULL;
};

#endif /* freesasaext_h */
