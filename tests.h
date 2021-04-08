#ifndef TESTS_INC
#define TESTS_INC

#include <vector>
#include <string>
#include <unordered_map>
#include <list>

#include "tools.h"
#include "index_reads.h"
#include "build_graph.h"

std::vector<std::vector<int>> true_adjMatrix(std::string tag, std::list<std::string> cloud);

void rapid_check();

#endif