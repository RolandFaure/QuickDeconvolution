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

void build_graph(short minCommonKmers, std::string tag, long int tagCloud, const std::vector<std::vector<long long int>> &readClouds, const std::vector <Read> &reads, const std::vector<std::set<long int>> &kmers, std::vector<int> &clusters);

void build_adj_matrix(short minCommonKmers, long int tagCloud, const std::vector <std::vector<long long int>> &readClouds, const std::vector <Read> &reads, const std::vector<std::set<long int>> &kmers, std::vector<std::vector<int>> &adjMatrix);

void fast_clustering(long int tagCloud, const std::vector <std::vector<long long int>> &readClouds, const std::vector <Read> &reads, const std::vector<std::set<long int>> &kmers, std::vector<int> &clusters);
#endif
