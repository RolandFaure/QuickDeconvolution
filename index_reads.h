#ifndef INDEX_READS_INCLUDED
#define INDEX_READS_INCLUDED

#include "tools.h"
#include "read.h"

#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>
#include <set>

void parse_reads(std::string fileReads, std::vector<std::vector<long long int>> &readClouds, std::vector <Read> &allreads, robin_hood::unordered_map <std::string, long int> &tagIDs, int min_length, int num_threads); //function to compile all reads before dispatching the work among the several threads

void compute_minimizers(int k, int h, int w, std::vector <Read> &allreads, int thread_id, int num_threads);

void index_kmers(int k, std::vector<std::vector<long int>> &kmers, std::vector <Read> &allreads, int thread_id, int num_threads);

void convert_kmers(std::vector<std::vector<long int>>& res, std::vector<std::vector<long int>>& input);

#endif //INDEX_READS_INCLUDED
