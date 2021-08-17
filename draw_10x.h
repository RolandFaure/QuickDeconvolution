#ifndef DRAW_10X_H_INCLUDED
#define DRAW_10X_H_INCLUDED

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <random>

//read_size      = 150       size of reads
//insert_size    = 400       # average size of insert
//gen_coverage   = 50        # genome sequencing coverage
//frag_coverage  = 0.2       # fragment sequencing coverage
//frag_size_min  = 20000     # minimum size of long fragments
//frag_size_max  = 70000     # maximum size of long fragments
//redundance     = 4         # average number of fragments sharing the same tag
void draw_fragments(int const read_size = 150, int insert_size = 400, int gen_coverage = 50, float frag_coverage = 0.2, int frag_size_min = 20000, int frag_size_max = 50000, float redundance = 4, float errorRate = 0, std::string fileGenome = "tests/genome.fasta", std::string fileReads = "tests/reads.fasta");

#endif // DRAW_10X_H_INCLUDED
