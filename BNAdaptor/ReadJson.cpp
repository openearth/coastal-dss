//
//  BNdemo.cpp
//  BNAdaptor
//
//  Created by Wiebke Jaeger on 25/09/2015.
//  Copyright (c) 2015 Wiebke Jaeger. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "json/json.h"

#include <boost/filesystem.hpp>
using namespace boost::filesystem;

Json::Value structure;
Json::Value BNStructure(void) {


    Json::Value file;
    if(exists("../workDir/BCnodes.json")){
        std::ifstream fileifs ("../workDir/info.json");
        fileifs>> file;
        std::cout << "info.json succesfully read" << std::endl;
    }
    else {
        std::cout << "ERROR: info.json is missing" << std::endl;
    }
    
    Json::Value BC;
    if(exists("../workDir/BCnodes.json")){
        std::ifstream BCifs ("../workDir/BCnodes.json");
        BCifs >> BC;
        std::cout << "BCnodes.json succesfully read" << std::endl;
    }
    else {
        std::cout << "ERROR: BCnodes.json is missing" << std::endl;
    }
        
    Json::Value R;
    if(exists("../workDir/Rnodes.json")){
        std::ifstream Rifs ("../workDir/Rnodes.json");
        Rifs >> R;
        std::cout << "Rnodes.json succesfully read" << std::endl;
    }
    else {
        std::cout << "ERROR: Rnodes.json is missing" << std::endl;
    }
    
    Json::Value H;
    if(exists("../workDir/Hnodes.json")){
        std::ifstream Hifs ("../workDir/Hnodes.json");
        Hifs >> H;
        std::cout << "Hnodes.json succesfully read" << std::endl;
    }
    else {
        std::cout << "ERROR: Hnodes.json is missing" << std::endl;
    }
    
    Json::Value C;
    if(exists("../workDir/Cnodes.json")){
        std::ifstream Cifs ("../workDir/Cnodes.json");
        Cifs >> C;
        std::cout << "Cnodes.json succesfully read" << std::endl;
    }
    else {
    std::cout << "WARNING: Cnodes.json is missing" << std::endl;
    }
    
	Json::Value M;
    if(exists("../workDir/Mnodes.json")){
        std::ifstream Mifs ("../workDir/Mnodes.json");
        Mifs >> M;
        std::cout << "Mnodes.json succesfully read" << std::endl;
    }
    else {
        std::cout << "WARNING: Mnodes.json is missing" << std::endl;
    }
    
	Json::Value Vul;
    if(exists("../workDir/vulnerabilities.json")){
            std::ifstream Vulifs ("../workDir/vulnerabilities.json");
        Vulifs >> Vul;
        std::cout << "vulnerabilities.json succesfully read" << std::endl;
    }
    else {
        std::cout << "WARNING: vulnerabilities.json is missing" << std::endl;
    }


    structure[0] = BC;
    structure[1] = R;
    structure[2] = H;
    structure[3] = C;
    structure[4] = M;
    structure[5] = file;
	structure[6] = Vul;
    return structure;
};