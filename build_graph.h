#ifndef BUILD_GRAPH_INCLUDED
#define BUILD_GRAPH_INCLUDED

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <fstream>
#include <chrono>
#include <list>

#include "tools.h"

std::vector<int> build_graph(short minCommonKmers, std::string tag, long int tagCloud, const std::vector<long long int>& readCloud, const std::vector <Read> &reads, std::vector<std::set<long int>> &kmers, std::vector<int> &clusters);

#endif
