//
//  testRandom.cpp
//  TPD_dummytarget
//
//  Created by Sebastian Swanson on 10/22/20.
//  Copyright © 2020 Sebastian Swanson. All rights reserved.
//

#include <cstdlib>
#include <iostream>

#include "msttypes.h"

int main() {
    MST::Structure();
    cout << "Test 'std::rand()" << endl;
    for (int i = 0; i < 10; i++) std::cout << std::rand() << std::endl;
    cout << "Test 'MstUtils::randInt()" << endl;
    for (int i = 0; i < 10; i++) std::cout << MstUtils::randInt(10) << std::endl;
    
    cout << "Now call srand(time(NULL))..." << endl;
    cout << "time: " << time(NULL) << endl;
    srand(time(NULL));
    cout << "Test 'std::rand()" << endl;
    for (int i = 0; i < 10; i++) std::cout << std::rand() << std::endl;
    cout << "Test 'MstUtils::randInt()" << endl;
    for (int i = 0; i < 10; i++) std::cout << MstUtils::randInt(10) << std::endl;
    return 0;
}
