#ifndef CLUSTER_GRAPH_H
#define CLUSTER_GRAPH_H

#include "tools.h"
#include <unordered_map>
#include <unordered_set>

void cluster_graph_chinese_whispers(robin_hood::unordered_map<long int, std::unordered_set<int>> &matching_tags, std::vector<int> &clusters, std::string &tag);

int find_threshold(float proportion, std::vector<int> &strengths_of_links);
#endif
