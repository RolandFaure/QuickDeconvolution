#ifndef INDEX_READS_INCLUDED
#define INDEX_READS_INCLUDED

#include "tools.h"

#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <list>

void index_reads(int k, int w, std::string fileReads, robin_hood::unordered_map<Sequence, std::array<std::vector<Hit>, 4>, Sequence::HashFunction> &index, std::vector<std::pair<std::string,std::vector<long long int>>> &readClouds, std::vector <Read> &allreads);

#endif //INDEX_READS_INCLUDED
