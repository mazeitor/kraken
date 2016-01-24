/*
 * Copyright 2015, Oriol Mazariegos Canellas <oriol.mazariegoscanellas@ndm.ox.ac.uk> University of Oxford and Derrick Wood <dwood@cs.jhu.edu>
 * 
 * This file is part of the Kraken taxonomic sequence classification system.
 *
 * Kraken is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Kraken is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Kraken.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file kraken-filter.cpp
 * @name kraken-filter
 * @author oriol mazariegos canellas
 * @author carlos del ojo elias and nicholas sanderson as a collaborators
 * @date 25 March 2015
 * @brief translated script from kraken-filter perl script by Derrick Wood using boost library in c++ 
 */

//compile: g++ -std=c++11 -fpermissive kraken-filter.cpp -o kraken-filter.o
//execute: time ../scripts/kraken-filter.o --threshold 0.15 --prefix 'path_to_kraken/' < 'path_to_file.kraken.txt' > file.txt

/**************************************
*************** includes **************
**************************************/
#include <sys/timeb.h>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <getopt.h>
#include <unordered_map>
#include <iomanip>
#include <omp.h>
#include <algorithm>

#define strtk_no_tr1_or_boost
#include "strtk.hpp"

#define version "v2015-09-01"

using namespace std;

/**************************************
*********** global variables ********** 
**************************************/
unordered_map<int,string> rank_map;
unordered_map<int,string> name_map;
unordered_map<int,vector<int>> child_lists;
unordered_map<int,int> clade_counts;
unordered_map<int,int> taxo_counts;
string prefix;
int show_zeros = 0;
int seq_count = 0;

/**************************************
************** functions **************
**************************************/

bool sort_vector (int i,int j) { 
	return (clade_counts[i]>clade_counts[j]); 
}

/**
 *  @brief      Check if file exist
 *  @param      string path of the file
 *  @return     bool true or false
*/
inline bool existsFile (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

/**
 *  @brief      Generic split string line
 *  @param      string String to split
 *  @param      char Delimiter to do the split
 *  @param      vector<string> Partial result
 *  @return     vector<string>
*/
vector<string> &split(const string &str, char delim, vector<std::string> &elems) {
    stringstream ss(str);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}
/**
 *  @brief      Generic recursive split string line
 *  @param      string String to split
 *  @param      char Delimiter to do the split
 *  @return     vector<string>
*/
vector<string> split(const string &str, char delim) {
    vector<string> elems;
    split(str, delim, elems);
    return elems;
}


/**
 *  @brief      Fast split string line
 *  @param      string String to split
 *  @param      char Delimiter to do the split
 *  @return     strtk::std_string::token_list_type
*/
strtk::std_string::token_list_type splitc7(string str, string delim)
{
	strtk::std_string::token_list_type token_list;
	strtk::split(delim,str,std::back_inserter(token_list));
	return token_list;
}

/**
 *  @brief      Convert a generic type number to string
 *  @param      Number generic type number
 *  @return     string a string number
*/
template <typename T> string NumberToString ( T Number ){
	ostringstream ss;
	ss << Number;
	return ss.str();
}

/**
 *  @brief     
*/
string rank_code(string rank_value)
{
	if(rank_value.compare("species") == 0) return "S";
	if(rank_value.compare("genus") == 0) return "G";
	if(rank_value.compare("family") == 0) return "F";
	if(rank_value.compare("order") == 0) return "O";
	if(rank_value.compare("class") == 0) return "C";
	if(rank_value.compare("phylum") == 0) return "P";
	if(rank_value.compare("kingdom") == 0) return "K";
	if(rank_value.compare("superkingdom") == 0) return "D";

	return "-";
}

/**
 *  @brief     
*/
void dfs_summation(int node)
{
    	//vector<int> children = child_lists.at(node); 
    	vector<int> children = child_lists[node]; 
	int aux = 0;

	int child = 0;
	int size = children.size();
	
	for(int i=0; i<size; i++){

		child = children[i];

		dfs_summation(child);

		aux = 0;

		auto search1 = clade_counts.find(child);
		if(search1 != clade_counts.end()) {
			aux = search1->second;
		}

		auto search2 = clade_counts.find(node);
		if(search2 != clade_counts.end()) {
			search2->second += aux;
		}else{
			clade_counts[node] = aux;
		}
		
	}

	children.clear();

}


/**
 *  @brief     
*/
void dfs_report(int node, int depth)
{
	float ccvalue = 0;
	ccvalue = (float)clade_counts[node];
	if(ccvalue == 0 && !show_zeros)	return;

	auto searchTC = taxo_counts.find(node);

	cout << " " << fixed << setprecision(2);
	cout << (ccvalue == 0 ? 0 : ccvalue * 100 / seq_count) << "\t";
	cout << (ccvalue == 0 ? 0 : (int)ccvalue) << "\t";
	cout << (searchTC == clade_counts.end() ? 0 : searchTC->second) << "\t";
	cout << rank_code(rank_map[node]) << "\t";

	cout << node << "\t";
	for(int i=0; i<depth; i++){
		cout << "  ";
	}
	cout << name_map[node] << endl;

	vector<int> children = child_lists[node]; 

	if(children.size() > 0){
		sort( children.begin(),  children.end(), sort_vector);
		for (vector<int>::iterator child =  children.begin(); child !=  children.end(); ++child) {
			dfs_report(*child, depth + 1 );
        	}
	}
	children.clear();
}

/**
 *  @brief	Load taxonomy names from kraken db      
 */
int load_taxonomy_names(){
	//check if file exist
	bool exist = existsFile((prefix+"/taxonomy/names.dmp").c_str());
	if(!exist){
		cout << "STDERR " << (prefix+"/taxonomy/names.dmp").c_str() << "must supply database name!" << endl; 
		return 1;
	}

	//open file
	ifstream NAMES((prefix+"/taxonomy/names.dmp").c_str());

	//processing file line by line
	string line;
	ios_base::sync_with_stdio(false); 
	while (getline(NAMES, line)) {
		istringstream iss(line);

		int node_id;
		string name;
		string type;
		string trash;

		//split only 0 and 1 and 3 position
		line.erase (std::remove(line.begin(), line.end(), '\t'), line.end());
		strtk::parse(line,"|",node_id,name,trash,type);

		if(type.compare("") == 0){
			type = trash;
			trash = "";
		}	

		if (type.compare("scientific name") == 0) {
			name_map[node_id] = name;
		}
	}
	NAMES.close();
	
	return 0;
}

/**
 *  @brief	Load taxonomy nodes from kraken db      
 */
int load_taxonomy_nodes(){
	//check if file exist
	bool exist = existsFile((prefix+"/taxonomy/nodes.dmp").c_str());
	if(!exist){
		cout << "STDERR " << (prefix+"/taxonomy/nodes.dmp").c_str() << "must supply database name!" << endl; 
		return 1;
	}

	//open file
	ifstream NODES((prefix+"/taxonomy/nodes.dmp").c_str());

	//processing file line by line
	string line;
	ios_base::sync_with_stdio(false); 
	while (getline(NODES, line)) {
		int node_id;
		int parent_id;
		string rank;
		string token;

		//split only 0 and 1 and 2 position
		line.erase (std::remove(line.begin(), line.end(), '\t'), line.end());
		strtk::parse(line,"|",node_id,parent_id,rank );

		if (node_id == 1) {
			parent_id = 0;
		}

		auto search = child_lists.find(parent_id);
		if(search != child_lists.end()) {
			search->second.push_back(node_id);
		}else{
			child_lists[parent_id].push_back(node_id);
		}

		rank_map[node_id] = rank;
	}
	NODES.close();

	return 0;
}

/**
 *  @brief	Load taxonomy from kraken db      
 */
int load_taxonomy(){
	load_taxonomy_names();
	load_taxonomy_nodes();
}

/**
 *  @brief	Process the filtering task for input files    
 */
int report(){
	taxo_counts[0] = 0;
	string fields;

	string field0;
	string field1;
	int field2;
	
	ios_base::sync_with_stdio(false); 
	while (getline(std::cin, fields)) {
		try{
			strtk::parse(fields,"\t",field0,field1,field2);
			auto search = taxo_counts.find(field2);
			if(search != taxo_counts.end()) {
				search->second++ ;
			}else{
				taxo_counts[field2] = 1;
			}
	    		seq_count++;
		}catch(const std::exception& ex){
		
		}
	}

	//int classified_count = seq_count - taxo_counts[0];
	int classified_count = 0;
	auto search = taxo_counts.find(0);
	if(search != taxo_counts.end()) {
		classified_count = seq_count - taxo_counts.at(0);
	}

	for (auto& it:taxo_counts) {
		//clade_counts.insert({it.first,it.second}); 
		clade_counts.insert({it}); 
	}

	dfs_summation(1);

	for (auto& it:name_map) {
		auto search = taxo_counts.find(it.first);
		if(search == taxo_counts.end()) {
			taxo_counts[it.first] = 0;
		}
	}

	double value = ((clade_counts[0] * 100) / seq_count);
	//printf("%6.2f \t", value );
	cout << " " << fixed << setprecision(2) << value << "\t";
	cout << clade_counts[0] << "\t";
	cout << taxo_counts[0] << "\t";
	cout << "U" << "\t";
	cout << "0" << "\t";
	cout << "unclassified" << endl;

	dfs_report( 1, 0 );

	return 0;
}

/****************************** HELP *****************************/
void print_help(){

	printf ("KRAKEN-REPORT application HELP - %s\n",version);

	printf ("db - path to database\n");

}

int main (int argc, char **argv)
{

	int c;
	int error;

	while (1) {
		static struct option long_options[] =
		{
			//These options don’t set a flag.
			//We distinguish them by their indices. 
			{"db",    required_argument, 	0, 	'd'},
			{"help",  no_argument, 		0, 	'h'},
			{0}
		};
		//getopt_long stores the option index here. 
		int option_index = 0;

		c = getopt_long (argc, argv, "d:h:",long_options, &option_index);

		//Detect the end of the options.
		if (c == -1) { 
			print_help();
			return 0;
			//break;
		}

		switch (c) {
			case 'd':
				prefix = optarg;
				break;
			case 'h':
			default:
				print_help();
				return 0;
		}
	}

	//loading kraken db
	error = load_taxonomy();
	if(error == 1) return 1;

	//filtering all stdin inputs
	error = report();
	if(error == 1) return 1;

	//deallocating memory
	rank_map.clear();
	name_map.clear();
	child_lists.clear();
	clade_counts.clear();
	taxo_counts.clear();
}

