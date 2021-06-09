#ifndef OUTPUT_H
#define OUTPUT_H

#include "tools.h"
#include "read.h"

#include <fstream>

void output(std::string inputFile, std::string outputFile, std::vector<Read>& reads);

#endif // OUTPUT_H
