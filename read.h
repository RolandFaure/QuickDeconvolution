#ifndef READ_H
#define READ_H

#include <vector>
#include <thread>
#include <mutex>

#include "tools.h"

class Read
{

public:
    Read(int num_threads);

    unsigned int barcode;
    short barcode_extension; //here a number to deconvolve the reads
    Sequence sequence;
    Sequence sequence_paired; //if the read is paired, here is the other end

    void new_mini(long int idx, int thread_id);
    std::vector<std::vector<long int>> &get_minis();
    void get_minis_seq(int k, std::vector<Sequence> &res, int thread_id);
    std::vector<std::vector<int>> get_pos_minis();

    void compute_minimisers(int k, int h, int w);

    void reset_minis_seq_thread (int thread); //a function to clear get_minis_seq when not needed anymore, to free memory

private :

    // each long int represents the index of the minimizer in the index
    std::vector<std::vector<long int>> minis; //one vector for each thread (in the simplest case, juste a vector)

    //each int is the position of a minimizer on the sequence
    std::vector<std::vector<int>> minis_seq;
    std::vector<std::vector<int>> minis_seq_rev;

    std::vector<std::vector<int>> minis_paired_seq;
    std::vector<std::vector<int>> minis_paired_seq_rev;

};

#endif // READ_H
