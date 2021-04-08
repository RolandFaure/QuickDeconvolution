#ifndef MEASURE_DEFINED
#define MEASURE_DEFINED

#include <vector>
#include <string>
#include <chrono>
#include <unordered_map>
#include <list>

#include "measure_speed.h"
#include "build_graph.h"
#include "index_reads.h"

float measure_graph_building_time(int k, int w, std::string readsFile);
void systematic_times(int k);

#endif