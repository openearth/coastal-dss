//
//  main.cpp
//  BNAdaptor
//
//  Created by Wiebke Jaeger on 21/09/2015.
//  Copyright (c) 2015 Wiebke Jaeger. All rights reserved.
//

#include <stdio.h>
#include "smile.h"
#include "math.h"
#include <iostream>
//#include "../json.hpp"
#include "json/json.h"

//using json = nlohmann::json;

Json::Value BNsetup;
std::vector<std::string> ncpath;

Json::Value BNStructure(void);
void CreateNetwork(Json::Value,std::vector<std::string>);
void FillTables(Json::Value,std::vector<std::string>);
std::vector<std::string> GetNCDirectories(void);
// We are reading 2D data, a 6 x 12 grid.

int main()
{
	BNsetup = BNStructure();
	ncpath = GetNCDirectories();
	CreateNetwork(BNsetup,ncpath);
	system("pause");
    return(DSL_OKAY);
};

