#ifndef CLUSTER_GRAPH_H
#define CLUSTER_GRAPH_H

#include "tools.h"
#include <unordered_map>
#include <unordered_set>

void cluster_graph_chinese_whispers(std::vector<std::vector<int>> &adjMatrix, std::vector<int> &clusters, std::string &tag);
#endif
