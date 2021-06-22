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

#include "read.h"
#include "tools.h"

void thread_deconvolve(short minCommonKmers, robin_hood::unordered_map <std::string, long int> &tagIDs, const std::vector <std::vector<long long int>> &readClouds, std::vector <Read> &reads, const std::vector<std::vector<std::vector<long int>>> &kmers, int thread_id, int num_thread, std::string folderOut);
void build_graph(short minCommonKmers, std::string tag, long int tagCloud, const std::vector<std::vector<long long int>> &readClouds, std::vector <Read> &reads, const std::vector<std::vector<std::vector<long int>>> &kmers, std::vector<int> &clusters, std::string folderOut);


void build_adj_matrix(short minCommonKmers, long int tagCloud, const std::vector <std::vector<long long int>> &readClouds, std::vector <Read> &reads, const std::vector<std::vector<std::vector<long int>>> &kmers, std::vector<std::vector<int>> &adjMatrix);

void fast_clustering(long int tagCloud, const std::vector <std::vector<long long int>> &readClouds, std::vector <Read> &reads, const std::vector<std::vector<std::vector<long int>>> &kmers, std::vector<int> &clusters);
void find_reps(long int tagCloud, const std::vector <std::vector<long long int>> &readClouds, std::vector <Read> &reads, const std::vector<std::vector<std::vector<long int>>> &kmers, std::vector<int> &clusterReps, robin_hood::unordered_map <long int, int> &alreadySeenTags);

#endif
