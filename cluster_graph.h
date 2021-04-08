#ifndef CLUSTER_GRAPH_H
#define CLUSTER_GRAPH_H

#include "tools.h"
#include <unordered_map>

void cluster_graph(std::unordered_map<long int, std::list<int>> &matching_tags, std::vector<int> &clusters);

//functions to find connected components
void find_connected_components(std::vector<std::vector<int>> &adjMatrix, std::vector<int> &clusters);
void extend_cluster_from_seed(int seed, std::vector<std::vector<int>> &adjMatrix, std::vector<int> &clusters);

int find_threshold(float proportion, std::vector<int> &strengths_of_links);
#endif