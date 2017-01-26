//
//  GetNCDirectories.cpp
//  BNAdaptor
//
//  Created by Wiebke Jaeger on 13/11/2015.
//  Copyright © 2015 Wiebke Jaeger. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <boost/filesystem.hpp>

using namespace boost::filesystem;

std::vector<std::string> NCPATH;
std::string PATH;

std::vector<std::string> GetNCDirectories()
{
    path p("../trainingData");
	
    for (auto i = directory_iterator(p); i != directory_iterator(); i++)
    {
        PATH = "../trainingData/";
		std::string hidden(i->path().filename().string());
        if (is_directory(p) && hidden.find(".") != 0){
            PATH.append(i->path().filename().string());
            PATH.append("/");
            NCPATH.push_back(PATH);  
            PATH.erase();
        }
    continue;
    }
    return NCPATH;
}