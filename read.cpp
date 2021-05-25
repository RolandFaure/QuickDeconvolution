#include "read.h"

using std::vector;

Read::Read(int num_threads)
{
    for (int t = 0 ; t < num_threads ; t++){
        minis.push_back({});
    }
    for (int t = 0 ; t < num_threads ; t++){
        minis_seq.push_back({});
    }
}

void Read::new_mini(long idx, int thread_id){

    //lock the mutex : only one thread can appen to minis at the same time
    //mutex.lock();

    minis[thread_id].push_back(idx);

    //mutex.unlock();
}

std::vector<std::vector<long int>>& Read::get_minis(){
    return minis;
}

std::vector<std::vector<Sequence>>& Read::get_minis_seq(){
    return minis_seq;
}

void Read::reset_minis_seq_thread (int thread){
    minis_seq[thread] = {};
}

void Read::compute_minimisers(int k, int h, int w){

    sequence.minimisers(h,k,w, minis_seq);
    sequence.reverse_complement().minimisers(h,k,w,minis_seq); //append to minis_seq also the minimisers of the reverse complement

    //now the sequence is useless, free some space by deleting it
    sequence = Sequence();
}
