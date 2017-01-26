//
//  main.cpp
//  BNAdaptor
//
//  Created by Wiebke Jaeger on 21/09/2015.
//  Copyright (c) 2015 Wiebke Jaeger. All rights reserved.

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "smile.h"
#include "netcdf.h"
#include <stdlib.h>
#include "json/json.h"
#include <boost/filesystem.hpp>
							

using namespace std;
using namespace boost::filesystem;
void GetStateNames(Json::Value ELEMENT, DSL_idArray &someNames);
void AddArcs(Json::Value PARENTS, Json::Value SETUP, vector<int> &BCPAR, vector<int> &RPAR, vector<int> &HPAR, vector<int> &CPAR, vector<int> &MPAR, vector<int> &MEffPAR, int &CHILD, DSL_network &theNET);

Json::Value setup;
string moduledatasetfilesDIR("../workDir/");

void CreateNetwork(Json::Value setup, vector<string> NCpath) {

	/* *** Preliminary Definitions *** */
	DSL_network theNet;
	DSL_intArray parcoord;
	DSL_idArray StateNames;
	DSL_idArray stateNames;
	vector<int> BC; BC.resize(setup[0].size());
	vector<int> R; R.resize(setup[1].size());
	vector<int> H; H.resize(setup[2].size());
	vector<int> C; C.resize(setup[3].size());
	vector<int> M; M.resize(setup[4].size());
	vector<int> MEff = M; 

	struct recinfo {  
		vector<int> gridid;
		vector<int> areaid;
		vector<int> receptorid;
	};
	vector<recinfo> allR;

	string BCfilenm = "hazardbc.nc";
	string Hfilenm = "hazard.nc";
	int bctime = setup[5]["hazardbc_time"].asInt(); 
	int bcstations = setup[5]["hazardbc_stations"].asInt(); 
	int htime = setup[5]["hazard_time"].asInt(); 
	int hstations;
	if(setup[5]["grid_type"] == "unstructured"){
		hstations = setup[5]["hazard_stations"].asInt(); 
		cout << "Grid recognized as unstructured. 'hazard_stations' = " << hstations << endl;
	} else if (setup[5]["grid_type"] == "structured") {
		hstations = setup[5]["hazard_row"].asInt()*setup[5]["hazard_col"].asInt();
		cout << "Grid recognized as structured. 'row' = "<< setup[5]["hazard_row"].asInt() << " and 'col' = " << setup[5]["hazard_col"].asInt()  << endl;
	} else {
		hstations = 0;
		cout << "ERROR: Undefined grid type." << endl;
	}

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

	/* *** Create Structure *** */

	// BOUNDARY CONDITIONS
	int bcind = 0;
	for (Json::ValueIterator itr =setup[0].begin(); itr != setup[0].end(); itr++){
		Json::Value element = *itr;
		// Add nodes
		string nm = element["title"].asString();
		const char* name = nm.c_str();
		if(element["nodetype"]!="deterministic"){
			BC[bcind] = theNet.AddNode(DSL_CPT,name);
		}
		else {
			BC[bcind] = theNet.AddNode(DSL_TRUTHTABLE,name);
		}
		// Add states
		GetStateNames(element,StateNames);
		theNet.GetNode(BC[bcind])->Info().Screen().color = 12759296;
		theNet.GetNode(BC[bcind])->Definition()->SetNumberOfOutcomes(StateNames);
		StateNames.Flush();
		// Add arcs (if parents exist)
		if (element["parents"] != "None"){
			AddArcs(element["parents"],setup,BC,R,H,C,M,MEff,BC[bcind],theNet);
		}
		bcind ++;
	}

   
	// MEASURES
	// Measures
	int mind = 0;
	for (Json::ValueIterator itr = setup[4].begin(); itr != setup[4].end(); itr++){
		Json::Value element = *itr;
		// Add node
		string nm = element["title"].asString();
		// cout << nm;
		const char* name = nm.c_str();
		M[mind] = theNet.AddNode(DSL_TRUTHTABLE,name);
		// Add states
		GetStateNames(element,StateNames);
		theNet.GetNode(M[mind])->Info().Screen().color = 47981;
		theNet.GetNode(M[mind])->Definition()->SetNumberOfOutcomes(StateNames);
		StateNames.Flush();
        
		// Add arcs (if parents exist)
		if (element["parents"] != "None"){
			AddArcs(element["parents"],setup,BC,R,H,C,M,MEff,M[mind],theNet);
		}
        
		// Add effectiveness node if effectiveness is < 1
		if(element["effectiveness"].asFloat() < 1){
            
			MEff[mind] = theNet.AddNode(DSL_CPT,"Implementation [%]");
			if (element["bins"][0] == "Yes") {
				StateNames.Add("Effective");
				StateNames.Add("Not_effective");
			}
			else if(element["bins"][0] == "No") {
				StateNames.Add("Not_effective");
				StateNames.Add("Effective");
			}
			else {
				cout<< "ERROR: Measure node does not exclusively have Yes and No as bins." << endl;
			}
			theNet.GetNode(MEff[mind])->Info().Screen().color = 12976103;
			theNet.GetNode(MEff[mind])->Definition()->SetNumberOfOutcomes(StateNames);
			StateNames.Flush();
			theNet.AddArc(M[mind],MEff[mind]);
		}
        
		mind++;
	}

	// RECEPTORS
	int rind = 0;
	for (Json::ValueIterator itr =setup[1].begin(); itr != setup[1].end(); itr++){
		Json::Value element = *itr;
		if(element["type"] != "latent"){
			// Add nodes
			string nm = element["title"].asString();
			const char* name = nm.c_str();
			R[rind] = theNet.AddNode(DSL_CPT,name);
			// Add states
			GetStateNames(element,StateNames);
			theNet.GetNode(R[rind])->Info().Screen().color = 11184810;
			theNet.GetNode(R[rind])->Definition()->SetNumberOfOutcomes(StateNames);
			StateNames.Flush();
			// Add arcs (if parents exist)
			if (element["parents"] != "None"){
				AddArcs(element["parents"],setup,BC,R,H,C,M,MEff,R[rind],theNet);
			}
		}
		rind++;
	}




	// HAZARDS
	int hind = 0;
	for (Json::ValueIterator itr =setup[2].begin(); itr != setup[2].end(); itr++){
		Json::Value element = *itr;
		// Add node
		string nm = element["title"].asString();
		const char* name = nm.c_str();
		H[hind] = theNet.AddNode(DSL_CPT,name);
		// Add states
		GetStateNames(element,StateNames);
		theNet.GetNode(H[hind])->Info().Screen().color = 4554202;
		theNet.GetNode(H[hind])->Definition()->SetNumberOfOutcomes(StateNames);
		StateNames.Flush();
		// Add arcs (if parents exist) 
		if (element["parents"] != "None"){
			AddArcs(element["parents"],setup,BC,R,H,C,M,MEff,H[hind],theNet);
		}
		hind++;
	}
    cout << "Here" << endl;

    

	// Impacts
	int cind = 0;
	for (Json::ValueIterator itr =setup[3].begin(); itr != setup[3].end(); itr++){
		Json::Value element = *itr;
		// Add node
		string nm = element["title"].asString();
		//                cout << nm;
		const char* name = nm.c_str();
		C[cind] = theNet.AddNode(DSL_CPT,name);
		// Add states
		GetStateNames(element,StateNames);
		theNet.GetNode(C[cind])->Info().Screen().color = 50687;
		theNet.GetNode(C[cind])->Definition()->SetNumberOfOutcomes(StateNames);
		StateNames.Flush();
		// Add arcs (if parents exist)
		if (element["parents"] != "None"){
			AddArcs(element["parents"],setup,BC,R,H,C,M,MEff,C[cind],theNet);
		}
		cind++;
	}


	/* *** Writing Tables *** */

    
    
    
	// Boundary Conditions
	bcind = 0;
	for (Json::ValueIterator itr =setup[0].begin(); itr != setup[0].end(); itr++){
		Json::Value element = *itr;
		// Add CPT
		if(element["nodetype"]!="deterministic"){
			int nobins = element["nobins"].asInt();
			int exp = 0;
			DSL_sysCoordinates theCoordinates(*theNet.GetNode(BC[bcind])->Definition());
			for(size_t i=0; i < NCpath.size(); i++) {
				string BCpath = NCpath[i];
				string filenm;
				for (auto i = directory_iterator(BCpath); i != directory_iterator(); i++)
				{
					filenm = i->path().filename().string();
					if(filenm != "hazard.nc")
					{
						BCpath.append(filenm);
					}
				}
				//cout << BCpath << endl;
				const char *BCPATH = BCpath.c_str();
#define BC_FILE_NAME BCPATH
				int ncid, varid;
				float *BC_data_in;
				BC_data_in = (float *) malloc(sizeof(float)*bctime*bcstations);
				/* Loop indexes, and error handling. */
				int retval;
				/* Open the file. NC_NOWRITE tells netCDF we want read-only access
				* to the file.*/
				if ((retval = nc_open(BC_FILE_NAME, NC_NOWRITE, &ncid)))
					ERR(retval);
				string ID = element["ID"].asString();
				const char *varname = ID.c_str();
				/* Get the varid of the data variable, based on its name. */
				if ((retval = nc_inq_varid(ncid, varname, &varid)))
					ERR(retval);
				/* Read the data. */
				if ((retval = nc_get_var_float(ncid, varid, &BC_data_in[0])))
					ERR(retval);
				/* Close the file, freeing all resources. */
				if ((retval = nc_close(ncid)))
					ERR(retval);
				float value = BC_data_in[0];
				for(int j =1; j < bctime; j++){ // find value of interest, since there are more than 1 timestep in the netcdf
					if (BC_data_in[j] < value){
						value = BC_data_in[j];
					}
				}	
				if(nobins == element["bins"].size()-1){
					float minbin = element["bins"][0].asFloat();
					float maxbin = element["bins"][nobins].asFloat();
					if (value < minbin || value >= maxbin) {
						cout << "Value " << value << " is out of bin range for node " << element["ID"] << endl;
					}
					for(int coord = 0; coord < nobins; coord++){
						theCoordinates.GoTo(coord);
						float binL = element["bins"][coord].asFloat();
						float binU = element["bins"][coord + 1].asFloat();
						if (value >= binL && value < binU) {
							theCoordinates.UncheckedValue() = (exp * theCoordinates.UncheckedValue() + 1) / (exp+1);
						} else {
							theCoordinates.UncheckedValue() = (exp * theCoordinates.UncheckedValue()) / (exp+1);
						}
					}
					exp += 1;
				}
				else if(nobins == element["bins"].size()){
					float minbin = element["bins"][0].asFloat();
					float maxbin = element["bins"][nobins-1].asFloat();
					if (value < minbin || value > maxbin) {
						cout << minbin << "-" << maxbin << endl;
						cout << "Value " << value << " is out of bin range for node " << element["ID"] << endl;
					}
					for(int coord = 0; coord < nobins; coord++){
						theCoordinates.GoTo(coord);
						float bin = element["bins"][coord].asFloat();
						if (value == bin) {
							theCoordinates.UncheckedValue() = (exp * theCoordinates.UncheckedValue() + 1) / (exp+1);
						} else {
							theCoordinates.UncheckedValue() = (exp * theCoordinates.UncheckedValue()) / (exp+1);
						}
					}
					exp += 1;
				}
				else{
					cout << "Misspecification in bins of BCnode " << element["ID"] << endl;
				}
			}
		}
		bcind ++;
		cout << "Table of " << element["ID"] << " filled" << endl;
	}
    

	// Measures
	mind=0;
	for (Json::ValueIterator itr = setup[4].begin(); itr != setup[4].end(); itr++){
		Json::Value element = *itr;
		if(element["effectiveness"].asFloat() < 1){
			DSL_sysCoordinates theCoordinates (*theNet.GetNode(MEff[mind])->Definition());
			if(element["bins"][0] == "Yes"){
				theCoordinates.UncheckedValue() = element["effectiveness"].asFloat();
				theCoordinates.Next();
				theCoordinates.UncheckedValue() = 1 - element["effectiveness"].asFloat();
				theCoordinates.Next();
				theCoordinates.UncheckedValue() = 0;
				theCoordinates.Next();
				theCoordinates.UncheckedValue() = 1;
			} else if (element["bins"][0] == "No")
			{
				theCoordinates.UncheckedValue() = 1;
				theCoordinates.Next();
				theCoordinates.UncheckedValue() = 0;
				theCoordinates.Next();
				theCoordinates.UncheckedValue() = 1 - element["effectiveness"].asFloat();
				theCoordinates.Next();
				theCoordinates.UncheckedValue() = element["effectiveness"].asFloat();
			}	else{
				cout << "WARNING: 'bins' should have values 'Yes' and 'No'" << endl;
			}
		}	
		mind++;
		cout << "Table of " << element["ID"] << " filled" << endl;
	}


	// Receptors 
	rind = 0;
	for (Json::ValueIterator itr =setup[1].begin(); itr != setup[1].end(); itr++){
		Json::Value element = *itr;
		DSL_intArray Rdim = 0;
		int nochildbins = element["nobins"].asInt();
		vector <double> probR;
		probR.resize(nochildbins);
		vector <int> areaid;
		vector <int> receptorid;
		vector <int> gridid;
		const char* Rfiletype = ".txt";
		string Rfilenm = element["ID"].asString();
		Rfilenm.append(Rfiletype);
		string Rpath = moduledatasetfilesDIR;
		Rpath.append(Rfilenm);
		Rpath =  Rpath.c_str();
		ifstream receptorfile (Rpath);
		if (receptorfile.is_open()) {
			while(! receptorfile.eof()){
				int areaidt, receptoridt, grididt;
				receptorfile >> areaidt >> receptoridt >> grididt;
				areaid.push_back(areaidt);
				receptorid.push_back(receptoridt);
				gridid.push_back(grididt);
			} // writes file content to areaid, gridid and receptorid
			allR.push_back(recinfo());
			allR[rind].gridid = gridid;
			allR[rind].areaid = areaid;
			allR[rind].receptorid = receptorid;
		} else {
			cout << "ERROR: Could not open text file for Rnode " << element["ID"] << "." << endl;
		}
		if(element["parents"] == "None"){
			for (int it = 0; it < nochildbins; it++) {
				probR[it] = 0;
			}
			int it = 0;	
			while (it < allR[rind].areaid.size()) {
				probR[allR[rind].areaid[it]-1]++;
				if (it + 1 < allR[rind].receptorid.size()){
					while (allR[rind].receptorid[it] == allR[rind].receptorid[it+1] && it+2 < allR[rind].receptorid.size()) {
						it++; 
					}
				}
				it++;
			}
			double sum = 0;
			for (int i=0; i < nochildbins; i++) {
				sum += probR[i];
			}
			for (int i=0; i < nochildbins; i++) {
				probR[i]=probR[i]/sum;
			}
			Rdim.Add(nochildbins);
			DSL_intArray cptcoord;
			DSL_Dmatrix coordinates = Rdim;
			DSL_sysCoordinates theCoordinates(*theNet.GetNode(R[rind])->Definition());
			if(element["type"] != "latent"){
				for (int childbin=0; childbin < nochildbins; childbin++) {
					cptcoord.Add(childbin);
					theCoordinates.GoTo(coordinates.CoordinatesToIndex(cptcoord));
					theCoordinates.UncheckedValue() = probR[childbin];
					cptcoord.Delete(cptcoord.GetSize()-1);
				}
			}
		}
		else if(element["parents"] != "None" && element["type"] != "latent"){ 	// initialize parent coordinates with array of zeros for CPT
			Json::Value parentcategory(Json::arrayValue);
			int noparents=0;
			if(element["parents"].size()>1){
				cout << "ERROR: Rnodes cannot be affected by more than one measure. Consider combining several measure nodes into one node with several bins. Probability tables will not be calculated for Rnode " << element["ID"] <<"." << endl;
				break;
			} 
			for (Json::ValueIterator itr2 =element["parents"].begin(); itr2 != element["parents"].end(); itr2++){			
				noparents +=1;
				Json::Value parent = *itr2;
				// Find in which category the node is to refind node handle
				for (int category = 0; category < 5; category++){
					int parind = 0;
					for (Json::ValueIterator itr3 = setup[category].begin(); itr3 != setup[category].end(); itr3++){
						Json::Value member = *itr3;
						if (member["ID"] == parent) {
							Rdim.Add(member["nobins"].asInt());
							switch(category){
							case 0:
								cout <<"BCnodes cannot be parents of Rnodes." << endl;
								break;
							case 1:
								cout <<"Rnodes cannot be parents of Rnodes." << endl;
								break;
							case 2:
								cout <<"Hnodes cannot be parents of Rnodes." << endl;
								break;
							case 3:
								cout <<"Cnodes cannot be parents of Rnodes." << endl;
								break;
							case 4:
								parentcategory.append("M");
								break;
							default:
								cout <<"Error in Assigning parent" << endl;
							}
						}
						parind++;
					}
				}
			}
			element["parentcat"] = parentcategory;
			// here begins the code for cpt with parents.
			Rdim.Add(nochildbins);
			DSL_intArray cptcoord;
			DSL_Dmatrix coordinates = Rdim;
			DSL_sysCoordinates theCoordinates(*theNet.GetNode(R[rind])->Definition());
			int j =0;
			for(int k=0; k<Rdim[j]; k++){ // Zeroth index refers to "No".
				cptcoord.Add(k);
				for (int it = 0; it < nochildbins; it++) {
					probR[it] = 0;
				}
				if(k==0){
					int it = 0;	
					while (it < allR[rind].areaid.size()) {
						probR[allR[rind].areaid[it]-1]++;
						if (it + 1 < allR[rind].receptorid.size()){
							while (allR[rind].receptorid[it] == allR[rind].receptorid[it+1] && it+2 < allR[rind].receptorid.size()) {
								it++; 
							}
						}
						it++;
					}
					double sum = 0;
					for (int i=0; i < nochildbins; i++) {
						sum += probR[i];
					}
					for (int i=0; i < nochildbins; i++) {
						probR[i]=probR[i]/sum;
					}
				}
				else{
					areaid.clear();
					receptorid.clear();
					Rfilenm = element["ID"].asString();
					string Mnm = element["parents"][j].asString();
					Rfilenm.append("_"); 
					Rfilenm.append(Mnm);
					Rfilenm.append("_");
					Rfilenm.append(static_cast<ostringstream*>( &(ostringstream() << k+1) )->str());
					Rfilenm.append(Rfiletype);
					Rpath = moduledatasetfilesDIR;
					Rpath.append(Rfilenm);
					Rpath =  Rpath.c_str();
					ifstream receptorfile (Rpath);
					if (receptorfile.is_open()) {
						while(! receptorfile.eof()){
							int areaidt, receptoridt, grididt;
							receptorfile >> areaidt >> receptoridt >> grididt;
							areaid.push_back(areaidt);
							receptorid.push_back(receptoridt);
						} 
					} else {
						cout << "ERROR: Could not open text file" << Rpath << "." << endl;
					}
					int it = 0;	
					while (it < areaid.size()) {
						probR[areaid[it]-1]++;
						if (it + 1 < receptorid.size()){
							while (receptorid[it] == receptorid[it+1] && it+2 < receptorid.size()) {
								it++; 
							}
						}
						it++;
					}
					double sum = 0;
					for (int i=0; i < nochildbins; i++) {
						sum += probR[i];
					}
					for (int i=0; i < nochildbins; i++) {
						probR[i]=probR[i]/sum;
					}
				}
				for (int childbin=0; childbin < nochildbins; childbin++) {
					cptcoord.Add(childbin);
					theCoordinates.GoTo(coordinates.CoordinatesToIndex(cptcoord));
					theCoordinates.UncheckedValue() = probR[childbin];
					cptcoord.Delete(cptcoord.GetSize()-1);
				}		
				cptcoord.Delete(cptcoord.GetSize()-1);
			}
		}
		else{
			cout << "ERROR: Parents of " << element["ID"] << " misspecified." << endl;
		}
		cout << "Table of " << element["ID"] << " filled" << endl;	
		rind++;
	}

	// Hazards
	hind = 0;
	for (Json::ValueIterator itr =setup[2].begin(); itr != setup[2].end(); itr++){
		Json::Value element = *itr;
		DSL_intArray Hdim = 0;
		DSL_intArray Hexpdim = 0;
		DSL_intArray cptcoord = 0;
		if (element["parents"] != "None"){
			// initialize parent coordinates with array of zeros for CPT
			Json::Value parentcategory(Json::arrayValue);
			for (Json::ValueIterator itr2 =element["parents"].begin(); itr2 != element["parents"].end(); itr2++){
				Json::Value parent = *itr2;
				cptcoord.Add(1);
				// Find in which category the node is to refind node handle
				for (int category = 0; category < 5; category++){
					int parind = 0;
					for (Json::ValueIterator itr3 = setup[category].begin(); itr3 != setup[category].end(); itr3++){
						Json::Value member = *itr3;
						if (member["ID"] == parent) {
							Hdim.Add(member["nobins"].asInt());
							Hexpdim.Add(member["nobins"].asInt());
							switch(category){
							case 0:
								parentcategory.append("BC");
								break;
							case 1:
								parentcategory.append("R");
								break;
							case 2:
								parentcategory.append("H");
								break;
							case 3:
								parentcategory.append("C");
								break;
							case 4:
								parentcategory.append("M");
								break;
							default:
								cout <<"Error in Assigning parent";
							}
						}
						parind++;
					}
				}
			}
			element["parentcat"] = parentcategory;
		}
		// Adding CPTs for individual H node. Loop over files, then parents. So that for each file, the CPT can be written immediately.
		int rec_pos_in_cptcoord = 99;
		DSL_Dmatrix  Hexp = Hexpdim;
		DSL_sysCoordinates expcoords(Hexp);
		for (int k = 0; k <= Hexp.GetSize()-1; k++) {  // initialize experience matrix with 0.
			expcoords.UncheckedValue() = 0;
			expcoords.Next();
		}
		int nochildbins = element["nobins"].asInt();
		Hdim.Add(nochildbins);
		DSL_Dmatrix coordinates = Hdim;
		DSL_sysCoordinates theCoordinates(*theNet.GetNode(H[hind])->Definition());
		for(int j=0; j < NCpath.size(); j++) {
			string BCpath = NCpath[j];
			string filenm;
			for (auto j = directory_iterator(BCpath); j != directory_iterator(); j++)
			{
				filenm = j->path().filename().string();
				if(filenm != "hazard.nc")
				{
					BCpath.append(filenm);
				}
			}
			const char *BCPATH = BCpath.c_str();
#define BC_FILE_NAME BCPATH
			for (int i=0; i<element["parents"].size(); i++) {
				if (element["parentcat"][i].asString() == "BC"){
					/* This will be the netCDF ID for the file and data variable. */
					int ncid, varid;
					float *BC_data_in;
					BC_data_in = (float *) malloc(sizeof(float)*bctime*bcstations);
					int retval;
					if ((retval = nc_open(BC_FILE_NAME, NC_NOWRITE, &ncid)))
						ERR(retval);
					string parID = element["parents"][i].asString();
					const char *varname = parID.c_str();
					if ((retval = nc_inq_varid(ncid, varname, &varid)))
						ERR(retval);
					if ((retval = nc_get_var_float(ncid, varid, &BC_data_in[0])))
						ERR(retval);
					if ((retval = nc_close(ncid)))
						ERR(retval);
					float value = BC_data_in[0];
					for(int j =1; j < bctime; j++){ // find value of interest, if there is more than 1 timestep in the netcdf
						if (BC_data_in[j] < value){
							value = BC_data_in[j];
						}
					}
					for (Json::ValueIterator itr4 =setup[0].begin(); itr4 != setup[0].end(); itr4++){
						Json::Value parent = *itr4;
						if(parent["ID"] == element["parents"][i]){
							int nobins = parent["nobins"].asInt();
							if(nobins == parent["bins"].size()-1){
								float minval = parent["bins"][0].asFloat();
								float maxval = parent["bins"][nobins].asFloat();
								if (value < minval || value >= maxval) {cout << "Value " << value << " is out of bin range for parent node " << parent["ID"] << endl;}
								else{
									for(int bin = 0; bin < nobins ; bin++){
										float binU = parent["bins"][bin + 1].asFloat();
										float binL = parent["bins"][bin].asFloat();
										if (value < binU && value >= binL) {cptcoord[i] = bin;}
									}
								}
							}
							else if(nobins == parent["bins"].size()){
								float minval = parent["bins"][0].asFloat();
								float maxval = parent["bins"][nobins-1].asFloat();
								if (value < minval || value > maxval) {cout << "Value " << value << " is out of bin range for parent node " << parent["ID"] << endl;}
								else{
									for(int bin = 0; bin < nobins ; bin++){
										if(value == parent["bins"][bin].asFloat()) {cptcoord[i] = bin;}
									}
								}
							}
						}
					}
				} else if (element["parentcat"][i] == "M"){
					string filenm_t;
					cptcoord[i] = 0;
					for (auto k = directory_iterator(NCpath[j]); k != directory_iterator(); k++)
					{
						filenm_t = k->path().filename().string();
						std::size_t extension = filenm_t.find(".");
						std::size_t first_ = filenm_t.find_first_of("_");
						std::size_t second_ = filenm_t.find_last_of("_");
						string measure = filenm_t.substr(first_+1,(second_-1)-first_);
						string bin = filenm_t.substr(second_+1,(extension-1)-second_);
						if (element["parents"][i] == measure)
							cptcoord[i] = stoi(bin)-1;
					}
					//Make if-loop that checks if there is an underscore in file name here and set coordinate to that value. Fix Hdim.

				}
				else if (element["parentcat"][i] == "R"){
					rec_pos_in_cptcoord = i;
				} else {
					cout << "ERROR: Parent nodes of " <<  element["ID"] << "misspecified." << endl;
				}
			}
			//        open hazard.nc in same directory and extract data
			string Hpath = NCpath[j];
			Hpath.append(Hfilenm);            
			string childIDstr = element["ID"].asString();
			istringstream childID (childIDstr);
			string hazardID;
			string receptorID;
			getline(childID, hazardID,'_');
			getline(childID, receptorID);
			const char *HPATH = Hpath.c_str();
#define H_FILE_NAME HPATH
			int ncid, varid;
			float *H_data_in;
			H_data_in = (float *) malloc(sizeof(float)*htime*hstations);
			int retval;
			if ((retval = nc_open(H_FILE_NAME, NC_NOWRITE, &ncid)))
				ERR(retval);
			const char *varname = hazardID.c_str();
			if ((retval = nc_inq_varid(ncid, varname, &varid)))
				ERR(retval);
			if ((retval = nc_get_var_float(ncid, varid, &H_data_in[0])))
				ERR(retval);
			if ((retval = nc_close(ncid)))
				ERR(retval);
			vector <float> hazval;
			for (int s = 0; s < hstations; s++) {
				hazval.push_back(H_data_in[s]);
			}
			/* compute probabilities & write cpt */
			vector <float> childbin_sums;
			vector<float> mean_hazval_at_rec;
			int rind = 0;
			for (Json::ValueIterator itr5 =setup[1].begin(); itr5 != setup[1].end(); itr5++){
				Json::Value parent = *itr5;
				if (parent["ID"] == receptorID) {
					int noareas = parent["nobins"].asInt();
					for (int area = 1; area <= noareas; area++) {				
						cptcoord[rec_pos_in_cptcoord] = area-1;
						int ncit = 0;
						int uniqueID = 0;
						mean_hazval_at_rec.clear();
						if(parent["bins"][area-1].asString() =="SafeZone"){
							mean_hazval_at_rec.push_back(0);
						}
						while (ncit < allR[rind].receptorid.size()){
							if (allR[rind].areaid[ncit] == area) {
								int dublicates = 1;
								mean_hazval_at_rec.push_back(hazval[allR[rind].gridid[ncit]-1]);
								ncit++;
								if(element["aggregation"] == "min"){
									if (ncit < allR[rind].receptorid.size()){
										while (allR[rind].receptorid[ncit-1] == allR[rind].receptorid[ncit] && ncit +1 < allR[rind].receptorid.size()) {
											if(mean_hazval_at_rec[uniqueID] > hazval[allR[rind].gridid[ncit]-1]){
												mean_hazval_at_rec[uniqueID] = hazval[allR[rind].gridid[ncit]-1];
											}
											ncit++;
										}
									}
								}
								else if(element["aggregation"] == "max"){
									if (ncit < allR[rind].receptorid.size()){
										while (allR[rind].receptorid[ncit-1] == allR[rind].receptorid[ncit] && ncit +1 < allR[rind].receptorid.size()) {
											if(mean_hazval_at_rec[uniqueID] < hazval[allR[rind].gridid[ncit]-1]){
												mean_hazval_at_rec[uniqueID] = hazval[allR[rind].gridid[ncit]-1];
											}
											ncit++;
										}
									}
								}
								else { //i.e. "aggregation" not specified or "mean"
									if (ncit < allR[rind].receptorid.size()){
										while (allR[rind].receptorid[ncit-1] == allR[rind].receptorid[ncit] && ncit +1 < allR[rind].receptorid.size()) {
											mean_hazval_at_rec[uniqueID] += hazval[allR[rind].gridid[ncit]-1];
											dublicates++;
											ncit++;
										}
									}
									mean_hazval_at_rec[uniqueID] /= dublicates;
								}
								uniqueID++;
							} else {ncit++;}
						}
						// compute fractions per bin here
						for (int childbin = 0; childbin < nochildbins; childbin++) {
							childbin_sums.push_back(0);
						}
						float minval = element["bins"][0].asFloat();
						int back  = element["nobins"].asInt();
						float maxval = element["bins"][back].asFloat();
						for (int pos = 0; pos < mean_hazval_at_rec.size(); pos++) {
							if (mean_hazval_at_rec[pos] < minval || mean_hazval_at_rec[pos] >= maxval){
								cout << "Value " << mean_hazval_at_rec[pos] << " is out of bin range for node " << element["ID"] << endl;
							} else {
								int childbin = 0;
								while (childbin < nochildbins){
									float binL = element["bins"][childbin].asFloat();
									float binU = element["bins"][childbin+1].asFloat();
									if (mean_hazval_at_rec[pos] >= binL && mean_hazval_at_rec[pos] < binU) { // taken out: cout << mean_hazval_at_rec[pos] << endl;
										childbin_sums[childbin] += 1;
										break;
									} else {childbin++;}
								}
							}
						}
						float sum_elem = 0;
						for (int childbin = 0; childbin < nochildbins; childbin++) { //this duplicate is wanted
							sum_elem += childbin_sums[childbin];
						}
						if(sum_elem != 0){
							float childbin_fractions;
							int exp = 0;
							for (int childbin = 0; childbin < nochildbins; childbin++) {
								childbin_fractions = childbin_sums[childbin]/sum_elem;
								expcoords.GoTo(Hexp.CoordinatesToIndex(cptcoord)); //coordinate of experience table
								cptcoord.Add(childbin);
								theCoordinates.GoTo(coordinates.CoordinatesToIndex(cptcoord)); // coordiante of cpt
								exp = expcoords.UncheckedValue();
								theCoordinates.UncheckedValue() = (theCoordinates.UncheckedValue()*exp +childbin_fractions) / (exp+1);
								cptcoord.Delete(cptcoord.GetSize()-1); // Deletes last entry again.
							}
							expcoords.UncheckedValue() += 1;
							childbin_sums.clear();
						}
					}
				}
				rind++;
			}
		}
		hind++;
		cout << "Table of " << element["ID"] << " filled" << endl;
	}

	// Impacts
	cind = 0;
	for (Json::ValueIterator itr =setup[3].begin(); itr != setup[3].end(); itr++){
		Json::Value element = *itr;
		if(element["type"].asString() == "categorical") {
			DSL_intArray Mdim;
			DSL_intArray MPardim;
			for (Json::ValueIterator itr2 =element["parents"].begin(); itr2 != element["parents"].end(); itr2++){
				Json::Value parent = *itr2;
				for (int category = 0; category < 5; category++){
					for (Json::ValueIterator itr3 = setup[category].begin(); itr3 != setup[category].end(); itr3++){
						Json::Value member = *itr3;
						if (member["ID"] == parent) {
							Mdim.Add(member["nobins"].asInt());
							MPardim.Add(member["nobins"].asInt());
						}
					}
				}
			}
			DSL_Dmatrix ParCoordinates(MPardim);
			Mdim.Add(element["nobins"].asInt()); // has now right number of dimensions and right dimensions
			DSL_Dmatrix coordinates = Mdim; // makes a matrix of system coordinates
			DSL_sysCoordinates theCoordinates(*theNet.GetNode(C[cind])->Definition());
			vector <int> vind; // makes the mapping a 1D array which can be transformed to coordinates using SMILE methods.
			for(int i =0; i < setup[6][element["ID"].asString()]["mapping"].size(); i++){
				vind.push_back(setup[6][element["ID"].asString()]["mapping"][i].asInt());
			}
			for(int i = 0; i < coordinates.GetSize(); i++){
				theCoordinates.UncheckedValue() = 0;	
				theCoordinates.Next();
			}	
			int items = MPardim.NumItems();
			for(int i = 0; i < vind.size(); i++){
				ParCoordinates.IndexToCoordinates(i,MPardim);
				DSL_intArray MappingCoord;
				for (int j = 0; j < items; j++){
					MappingCoord.Add(MPardim[j]);
				}
				MappingCoord.Add(vind[i]-1); // should be zero based
				theCoordinates.GoTo(coordinates.CoordinatesToIndex(MappingCoord));
				theCoordinates.UncheckedValue() = 1;
			}
		}
		cind++;
		cout << "Table of " << element["ID"] << " filled" << endl;
	}

	cout << "Process of filling tables completed" << endl;


	/* *** Save to file  *** */
	theNet.SimpleGraphLayout();
	const char* filename;
	const char* filetype = ".dsl";
	string filenm = "../output/";
	filenm.append(setup[5]["name"].asString());
	filenm.append(filetype);
	filename =  filenm.c_str();
	theNet.WriteFile(filename);

	cout << "BN saved as: " << filename << endl;
	};