#include "json/json.h"
#include "smile.h"
#include <iostream>



void GetStateNames(Json::Value ELEMENT, DSL_idArray &someNames)
{	
	const char* stname;
    int nobins = ELEMENT["nobins"].asInt();
    if (nobins == ELEMENT["bins"].size() - 1 && ELEMENT["bins"][1].isNumeric()) {
        for (int s=0; s < ELEMENT["bins"].size()-1; s++) {
            double stL = ELEMENT["bins"][s].asDouble();
            char stLchar[20];
            sprintf (stLchar,"%g", stL);
            std::string stLstr(stLchar);
            if (stLstr[0] == '-') {
                stLstr[0] = 'm';
            } else {
                stLstr.insert(0, "p");
            }
            double stU = ELEMENT["bins"][s+1].asDouble();
            char stUchar[20];
            sprintf (stUchar,"%g", stU);
            std::string stUstr(stUchar);
            if (stUstr[0] == '-') {
                stUstr[0] = 'm';
            } else {
                stUstr.insert(0, "p");
            }
            std::string st = stLstr  + "__" + stUstr;
            for (int i=0; i< st.size(); i++){
                if (st[i] == '.') {
                    st[i] = '_';
                }
				if (st[i] == '-') {
                    st[i] = 'm';
                }
				if (st[i] == '+') {
                    st[i] = 'p';
                }
            }
            stname = st.c_str();
            someNames.Add(stname);
        }
    }
    else if (nobins == ELEMENT["bins"].size() && ELEMENT["bins"][1].isNumeric()){
        for (int s=0; s < ELEMENT["bins"].size(); s++) {
            double stL = ELEMENT["bins"][s].asDouble();
            char stLchar[20];
            sprintf (stLchar,"%g", stL);
            std::string stLstr(stLchar);
            if (stLstr[0] == '-') {
                stLstr[0] = 'm';
            } else {
                stLstr.insert(0, "p");
            }
			std::string st = stLstr;
            for (int i=0; i< st.size(); i++){
                if (st[i] == '.') {
                    st[i] = '_';
                }
				if (st[i] == '-') {
                    st[i] = 'm';
                }
				if (st[i] == '+') {
                    st[i] = 'p';
                }
            }
            stname = st.c_str();
            someNames.Add(stname);
        }
    }
    else if (nobins == ELEMENT["bins"].size() && ELEMENT["bins"][1].isString()){
                for (Json::ValueIterator itr =ELEMENT["bins"].begin(); itr != ELEMENT["bins"].end(); itr++){
                    Json::Value state = *itr;
                    std::string st = state.asString();
                    stname = st.c_str();
					        for (int i=0; i< st.size(); i++){
                if (st[i] == ' ') {
                    st[i] = '_';
                }
            }
                    someNames.Add(stname);
                }
            }
    else {
        std::cout << "nobins and bins do not match for node" << ELEMENT["title"] << std::endl;
    }

}
