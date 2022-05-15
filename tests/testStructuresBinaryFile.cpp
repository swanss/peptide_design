//
//  testStructuresBinaryFile.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 11/27/21.
//  Copyright © 2021 Sebastian Swanson. All rights reserved.
//

#include <stdio.h>

//mst dependencies
#include "mstsystem.h"
#include "mstoptions.h"

//tpd dependencies
#include "utilities.h"
#include "coverage.h"
#include "benchmarkutilities.h"
//#include "termextension.h"
#include "secondarystructure.h"

int main(int argc, char *argv[]) {
    MstOptions op;
    op.setTitle("");
    op.addOption("bin_path", "path to a seed binary file.",true);
//    op.addOption("structure", "the name of a structure which will be loaded with detail",false);

    op.setOptions(argc, argv);
    
    // Variables provided by user
    string extfrag_bin_path = op.getString("bin_path");
//    string structureName = op.getString("structure","");
    
    // Read seeds file
    bool read = true;
    StructuresBinaryFile* extfrag_bin = new StructuresBinaryFile(extfrag_bin_path,read);
    extfrag_bin->scanFilePositions();
    
    return 0;
}
    
