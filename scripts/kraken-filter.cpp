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
 * @brief translated script from kraken-filter perl script by Derrick Wood applying multythreading and boost library in c++
 */

//compile: g++-4.7 -std=c++11 -Ofast kraken-filter-pthread.cpp -lpthread -o kraken-filter.o
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
#include <stdio.h>

#define strtk_no_tr1_or_boost
#include "strtk.hpp"
#include <pthread.h>

#define version "v2015-09-01"

using namespace std;



/**************************************
************** structures ************* 
**************************************/
struct thread_data{
	int tid;
	char *field;
};

/**************************************
*********** global variables ********** 
**************************************/
unordered_map<int,int> parent_map;
string prefix;
string input_file;
string output_file;
float threshold;
int MAX_THREADS = 16;
pthread_mutex_t count_mutex = PTHREAD_MUTEX_INITIALIZER;


/**************************************
************** functions **************
**************************************/

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
 *  @brief      Fast split string line
 *  @param      string String to split
 *  @param      char Delimiter to do the split
 *  @return     strtk::std_string::token_list_type
*/
strtk::std_string::token_list_type splitc7(string str, char delim)
{
	strtk::std_string::token_list_type token_list;
	strtk::split(delim,str,std::back_inserter(token_list));
	return token_list;
}

/**
 *  @brief	Load taxonomy from kraken db      
 */
int load_taxonomy(){

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
	while (getline(NODES, line)) {
		istringstream iss(line);

		int node_id;
		int parent_id;
		string token;

		//split only 0 and 1 position
		iss >> node_id >> token >> parent_id;

		if (node_id == 1) {
			parent_id = 0;
		}

		parent_map[node_id] = parent_id;
	}
	NODES.close();
	
	return 0;
}

/**
 *  @brief	Process the filtering task for input files    
 */
void *filtering(void *t){

	try{
		char const* delims = "\t";
		vector<string> splited_line;
		string code;
		string seqid;
		int called_taxon;
		string len;
		string hit_list;
		string taxid ;
		int ct;
		vector<string> hits; 

		unordered_map<string,int> hit_counts;
		unordered_map<int,int> hit_sums;
	
		string fields;

		ios_base::sync_with_stdio(false); 

		bool exit = false;
		while (true) {
			pthread_mutex_lock( &count_mutex );
			if(!getline(std::cin, fields)) {
				exit = true;
			}
			pthread_mutex_unlock( &count_mutex );


			if(exit == true){
				splited_line.clear();
				hits.clear(); 

				hit_counts.clear();
				hit_sums.clear();

				pthread_exit(NULL);
				return 0;
			}
		

			hit_counts.clear();

			strtk::parse(fields,"\t",code,seqid,called_taxon,len,hit_list);

			strtk::std_string::token_list_type token_list = splitc7(hit_list, ' ');
			strtk::std_string::token_list_type::iterator itr = token_list.begin();
			while (token_list.end() != itr){
				strtk::parse(string(itr->first, itr->second),":",taxid,ct);

				auto search = hit_counts.find(taxid);
				if(search != hit_counts.end()) {
					search->second += ct;
				}
				else {
					hit_counts[taxid] += ct;
				}

				++itr;
			}	

			hit_sums.clear();
			int total_kmers = 0;
			int total_unambig = 0;
			int total_hits = 0;

			int count;

			for (auto& it:hit_counts) {
				count = it.second;
	
				total_kmers += count;

				if (it.first.compare("A") != 0) {
					total_unambig += count;

					if (atoi(it.first.c_str()) > 0) {
						total_hits += count;

						int taxidAux = atoi(it.first.c_str());
						while (taxidAux > 0){
							auto search = hit_sums.find(taxidAux);
							if(search != hit_sums.end()) {
								search->second += count;
							}
							else {
								hit_sums[taxidAux] += count;
							}

							taxidAux = parent_map.at(taxidAux); 
						}				
					}
				}
			}

			float pct = 0;

			int new_taxon = called_taxon;

			while ( new_taxon > 0 )	{
				auto search = hit_sums.find(new_taxon);
				if(search != hit_sums.end()) {
					pct = (float)(search->second) / total_unambig;
				}
				if ( pct >= threshold - (1/100000) ){    
					break;
				}
				new_taxon = parent_map.at(new_taxon);
			}


			pthread_mutex_lock( &count_mutex );

			cout << (new_taxon > 0 ? "C" : "U") << "\t" << seqid << "\t" << new_taxon << "\t" << len << "\t" << "P=" << fixed << setprecision(3) << pct << "\t" << hit_list << endl;

			pthread_mutex_unlock( &count_mutex );

		}

		splited_line.clear();
		hits.clear(); 

		hit_counts.clear();
		hit_sums.clear();

		pthread_exit(NULL);
	}catch(const std::exception& ex){
		
	}
}

/****************************** HELP *****************************/
void print_help(){

	printf ("KRAKEN-FILTER application HELP - %s\n",version);

	printf ("db - path to database\n");
	printf ("threshold - threshold value\n");
	printf ("threads - threads to execute in parallel\n");

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
			{"db",    	required_argument, 	0, 	'd'},
			{"threshold",   required_argument, 	0, 	'r'},
			{"threads",   	required_argument, 	0,	't'},
			{"help", 	no_argument, 		0, 	'h'},
			{	0,			0,	0,	  0}
		};
		//getopt_long stores the option index here. 
		int option_index = 0;

		c = getopt_long (argc, argv, "d:r:h:t:",long_options, &option_index);

		//Detect the end of the options.
		if (c == -1) { 
			print_help();
			return 0;
		}

		switch (c) {
			case 'd':
				prefix = optarg;
				break;
			case 'r':
				threshold = stof(optarg);
				break;
			case 't':
				MAX_THREADS = atoi(optarg);
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

	//instantiate threads array
	pthread_t *thread_id = new pthread_t[MAX_THREADS];

	//create threads
	for(int i = 0; i < MAX_THREADS; i++)
	{
		int rc = pthread_create(&thread_id[i], NULL, filtering, NULL);
		if (rc != 0)
		{
			cerr << "Error creating thread tid=" << rc << endl;
		}
	}

	//joint each threat to main
	for(int i = 0; i < MAX_THREADS; i++)
	{
		pthread_join(thread_id[i], NULL);
	}

	//deallocating structure thread
	delete [] thread_id;

	if(error == 1) return 1;
}
