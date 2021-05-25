#include "index_reads.h"

#include <chrono>

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::pair;
using std::string;
using std::ofstream;
using std::ifstream;
using std::array;
using std::set;
using robin_hood::unordered_map;
using namespace std::chrono;



//to compress sequencing data, reads and kmers are represented by vectors of bool (encapsulated in Sequence)

//function parsing all the reads present in the file
void parse_reads(std::string fileReads, std::vector<std::vector<long long int>> &readClouds, std::vector <Read> &allreads, unordered_map <string, long int> &tagIDs, int num_threads){

    auto t0 = high_resolution_clock::now();
    char format = '@'; //a character to keep track of whether the input file is a fasta or a fastq
    if (fileReads.substr(fileReads.size()-6,6) == ".fasta"){
        format = '>';
    }

    ifstream in(fileReads);
    if (!in){cout << "problem reading files in index_reads, while trying to read " << fileReads << endl;}

    string nameofsequence;
    long long int sequenceID = 0; //instead of storing the name of sequence in a string we're going to store them in long long int, it will be more efficient

    string tag;
    bool next = false;

    long int tagID;

    string line;
    bool notag = false; //bool alerting when there is no tag attached to a read
    while(getline(in, line)){

        if (line[0] == format){

            Read r(num_threads);

            //here looking at the name of sequence and the tag
            nameofsequence = line.erase(0,1);
            tag = get_tag(nameofsequence, format); //this tag is a string, as contained in a fasta/q: we're going to convert it into a long int, this will be much more efficient

            if (tag == ""){}
            else if (tagIDs.find(tag) == tagIDs.end()){ // if true, tag is not in the keys of tagIDs

                tagIDs[tag] = readClouds.size(); //from now on this tag ID will be associated to this tag
                tagID = readClouds.size();
                readClouds.push_back({sequenceID});

                r.barcode = tagID;
                allreads.push_back(r);

                sequenceID ++;
                next = true;
            }
            else{
                tagID = tagIDs[tag];
                readClouds[tagID].push_back(sequenceID);
                r.barcode = tagID;
                allreads.push_back(r);

                sequenceID ++;
                next = true;
            }

        }
        else if (next) {
            //here looking at the sequence itself
            next = false;

            allreads[allreads.size()-1].sequence = Sequence(line);
        }
    }
    auto t1 = high_resolution_clock::now();
    cout << "Finished parsing all the reads in " << duration_cast<seconds>(t1-t0).count() << "s, for a total of " << sequenceID << " reads" << endl;
}

//function computing minimisers (one thread)
void compute_minimizers(int k, int h, int w, std::vector <Read> &allreads, int thread_id, int num_threads){

    for (int r = thread_id ; r < allreads.size() ; r+=num_threads){
        allreads[r].compute_minimisers(k, h, w);

        if ((r-thread_id)%(100000*num_threads) == 0){
            cout << "Thread " << thread_id << " computed " << double(r)/double(allreads.size())*100 << "% of its minimizers" << endl;
        }
    }

}

//this function partially fills the index (kmers) with one thread
void index_kmers(vector<vector<long int>> &kmers, vector <Read> &allreads, int thread_id, int num_threads){ //kmers is a vector of vector (and not vec of vec of vec) becuse we're focusing on 1 thread

    double total_mini_time = 0;
    double total_read_time = 0;
    double expansion_kmers_time = 0;
    long int count = 0;
    auto t0 = high_resolution_clock::now();

    //maps all minimizing kmers to their index in kmers
    unordered_map<Sequence,long long int, Sequence::HashFunction> index;

    for (long int r = 0 ; r < allreads.size() ; r++){

         //here we split the work among all thread: this thread just takes care of a few minimizing kmers


        for (Sequence mini : allreads[r].get_minis_seq()[thread_id]){


            auto tt0 = high_resolution_clock::now();
            if (index.find(mini) != index.end()){
                auto ttt0 = high_resolution_clock::now();
                long int place = index[mini];
                allreads[r].new_mini(place, thread_id);

                kmers[place].push_back(allreads[r].barcode); //each kmer contains the set of all barcodes containing this kmer
                auto ttt1 = high_resolution_clock::now();
                expansion_kmers_time += duration_cast<nanoseconds>(ttt1 - ttt0).count();
            }
            else {

                auto ttt0 = high_resolution_clock::now();
                index[mini] = kmers.size();
                allreads[r].new_mini(kmers.size(), thread_id);

                kmers.push_back({allreads[r].barcode});
                auto ttt1 = high_resolution_clock::now();
                expansion_kmers_time += duration_cast<nanoseconds>(ttt1 - ttt0).count();
            }

            auto tt1 = high_resolution_clock::now();
            total_read_time += duration_cast<nanoseconds>(tt1 - tt0).count();

        }

        count ++;
        if (count % 50000 == 0){
            cout << "Thread " << thread_id << " processed " << double(count)/double(allreads.size())*100 << "% of all sequences" << endl;
        }

    }

    auto t1 = high_resolution_clock::now();

    cout << "In thread " << thread_id << " of index_reads, handling the memory of kmers took me reads " << total_read_time/1000000000 << "s out of " << duration_cast<microseconds>(t1 - t0).count()/1000000 << "s in total" <<  endl;
}

//function returning all minimizers of a sequence and their positions knowing k and the windowsize


