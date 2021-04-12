#ifndef INDEX_READS_INCLUDED
#define INDEX_READS_INCLUDED

#include "tools.h"

#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>

void index_reads(int k, int w, std::string fileReads, robin_hood::unordered_map<Sequence, std::vector<std::vector<Hit>>, Sequence::HashFunction> &index, std::vector<std::vector<long long int>> &readClouds, std::vector <Read> &allreads);

#endif //INDEX_READS_INCLUDED
