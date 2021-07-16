#ifndef MEASURE_DEFINED
#define MEASURE_DEFINED

#include <vector>
#include <string>
#include <chrono>
#include <unordered_map>
#include <list>

#include "build_graph.h"
#include "index_reads.h"
#include "output.h"

float measure_graph_building_time(int k, int h, int w, int num_threads, std::string readsFile, std::string folderOut, std::string fileOut);

void systematic_times(int k);

#endif
