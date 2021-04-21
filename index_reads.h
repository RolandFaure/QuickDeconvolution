#ifndef INDEX_READS_INCLUDED
#define INDEX_READS_INCLUDED

#include "tools.h"

#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>
#include <set>

void index_reads(int k, int h, int w, std::string fileReads, std::vector<std::vector<long int>> &kmers, std::vector<std::vector<long long int>> &readClouds, std::vector <Read> &allreads);

#endif //INDEX_READS_INCLUDED
