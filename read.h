#ifndef READ_H
#define READ_H

#include <vector>
#include <thread>
#include <mutex>
#include <shared_mutex>

#include "tools.h"

class Read
{

public:
    Read(int num_threads);

    long int barcode;
    Sequence sequence;

    void new_mini(long int idx, int thread_id);
    std::vector<std::vector<long int>> &get_minis();
    std::vector<std::vector<Sequence>>& get_minis_seq();

    void compute_minimisers(int k, int h, int w);

    void reset_minis_seq_thread (int thread); //a function to clear get_minis_seq when not needed anymore, to free memory

private :

    // each long int represents the index of the minimizer in the index
    std::vector<std::vector<long int>> minis; //one vector for each thread (in the simplest case, juste a vector)

    //each Sequence is a minimizer
    std::vector<std::vector<Sequence>> minis_seq;

};

#endif // READ_H
