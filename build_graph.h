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
#include <boost/graph/adjacency_list.hpp>

#include "tools.h"

std::vector<int> build_graph(int k, int w, long int tagCloud, std::vector<long long int>& readCloud, std::vector <Read> &reads, std::unordered_map<Sequence, std::vector<std::vector<Hit>>, Sequence::HashFunction> &index);
#endif