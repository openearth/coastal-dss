#include "json/json.h"
#include "smile.h"
#include <iostream>
#include <vector>

using namespace std;

void AddArcs(Json::Value PARENTS, Json::Value SETUP, vector<int> &BCPAR, vector<int> &RPAR, vector<int> &HPAR, vector<int> &CPAR, vector<int> &MPAR, vector<int> &MEffPAR, int &CHILD, DSL_network &theNET){

	for (Json::ValueIterator itr2 =PARENTS.begin(); itr2 != PARENTS.end(); itr2++){
		Json::Value parent = *itr2;
		for (int category = 0; category < 5; category++){
			int parind = 0;
			for (Json::ValueIterator itr3 =SETUP[category].begin(); itr3 != SETUP[category].end(); itr3++){
				Json::Value member = *itr3;
				if (member["ID"] == parent) {
					//cout << member["ID"];
					switch(category){
					case 0:
						theNET.AddArc(BCPAR[parind],CHILD);
						break;
					case 1:
						theNET.AddArc(RPAR[parind],CHILD);
						break;
					case 2:
						theNET.AddArc(HPAR[parind],CHILD);
						break;
					case 3:
						theNET.AddArc(CPAR[parind],CHILD);
						break;
					case 4:
						if(member["effectiveness"].asFloat()==1){
							theNET.AddArc(MPAR[parind],CHILD);
						}
						else if(member["effectiveness"].asFloat() < 1 && member["effectiveness"].asFloat() > 0){
							theNET.AddArc(MEffPAR[parind],CHILD);
						}
						else {
							cout << "ERROR: Misspecification of effectiveness" << endl;
						}
						break;
					default:
						cout <<"Error in Assigning parent";
					}
				}
				parind++;
			}
		}
	}
}
