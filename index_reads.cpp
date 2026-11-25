#include "index_reads.h"

#include <chrono>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <stdexcept>

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
void parse_reads(
    std::string fileReads,
    std::vector<std::vector<long long int>>& readClouds,
    std::vector<Read>& allreads,
    robin_hood::unordered_map<std::string, long int>& tagIDs,
    int min_length,
    int num_threads,
    bool verbose)
{
    const auto started = std::chrono::steady_clock::now();
    const char format = sequence_format(fileReads);

    std::ifstream input_file;
    boost::iostreams::filtering_istream in;
    open_sequence_input(fileReads, input_file, in);

    long long int sequence_id = 0;
    long long int skipped_reads = 0;
    std::string previous_key;
    std::string tag;
    std::string line;
    std::vector<std::string> buffer;

    const auto parse_record = [&]() {
        if (buffer.empty()) {
            return;
        }
        if (!to_deconvolve(buffer, format, min_length, tag)) {
            ++skipped_reads;
            return;
        }

        const auto key = mate_key(buffer[0], format);
        if (sequence_id > 0 && key == previous_key) {
            allreads.back().sequence_paired = buffer[1];
            return;
        }

        const auto existing_tag = tagIDs.find(tag);
        long int tag_id;
        if (existing_tag == tagIDs.end()) {
            tag_id = static_cast<long int>(readClouds.size());
            tagIDs[tag] = tag_id;
            readClouds.push_back({sequence_id});
        } else {
            tag_id = existing_tag->second;
            readClouds[tag_id].push_back(sequence_id);
        }

        Read read(num_threads);
        read.sequence = buffer[1];
        read.barcode = tag_id;
        allreads.push_back(read);

        previous_key = key;
        ++sequence_id;
    };

    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == format) {
            parse_record();
            buffer = {line};
        } else {
            buffer.push_back(line);
        }
    }
    parse_record();

    const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::steady_clock::now() - started).count();
    if(verbose){
        std::cout << "Finished parsing in " << elapsed << "s. "
                  << sequence_id << " read entries to deconvolve, "
                  << skipped_reads << " records skipped.\n";
    }
}

//function computing minimisers of some reads(one thread)
void compute_minimizers(int k, int h, int w, std::vector <Read> &allreads, int thread_id, int num_threads, bool verbose){

    auto t0 = high_resolution_clock::now();

    for (int r = thread_id ; r < allreads.size() ; r+=num_threads){
        Read &re = allreads[r];
        re.compute_minimisers(k, h, w);

    //    if (r == 2989710){
    //        cout << "computed minimizers of " << re.sequence.str() <<" (barcode " << re.barcode << ") : " << endl;
    //        for (auto minit : re.get_pos_minis()){
    //            for (long int n : minit){
    //                cout << " - " << re.sequence.subseq(n,k).str() << " (pos minis : " << n << ")" << endl;
    //            }
    //        }
    //    }

        if ((r-thread_id)%(100000*num_threads) == 0 && verbose){
            cout << "Thread " << thread_id << " computed " << double(r)/double(allreads.size())*100 << "% of its minimizers" << endl;
        }
    }

    auto t1 = high_resolution_clock::now();

    if(verbose){
        cout << "Thread " << thread_id << " computed its minimizers in " << duration_cast<seconds>(t1-t0).count()<< "s" << endl;
    }
}

//this function partially fills the index (kmers) with one thread
void index_kmers(int k, vector<vector<long int>> &kmers, vector <Read> &allreads, int thread_id, int num_threads, bool verbose){ //kmers is a vector of vector (and not vec of vec of vec) becuse we're focusing on 1 thread

    std::mutex mutex;

    double total_read_time = 0;
    double expansion_kmers_time = 0;
    long long int count = 0;
    auto t0 = high_resolution_clock::now();

    int dividedsize = allreads.size()/num_threads;
    int asize = allreads.size();

    //maps all minimizing kmers to their index in kmers
    unordered_map<Sequence,long long int, Sequence::HashFunction> index;

    for (long int rr = 0 ; rr < asize ; rr++){

            int r = (rr + dividedsize*thread_id)%asize; //that is a line to encourage each thread to work on separate reads

            //here we split the work among all thread: this thread just takes care of a few minimizing kmers

            vector<Sequence> minis;
            Read &re = allreads[r]; //we will handle a reference to allreads[r] to avoid accessing allreads too much

            re.get_minis_seq(k, minis, thread_id); //this function only returns the kmers appropriate for this thread, they can all be appended

    //        if (thread_id == 0 && r == 95353){
    //            cout << "read 95353 (" << re.sequence.str() << "):";
    //            for (Sequence mini : minis) cout << mini.str() << ",";
    //            cout << endl;
    //        }
    //        if (thread_id == 1 && r == 95379){
    //            cout << "read 95352 (" << re.sequence.str() << "):";
    //            for (Sequence mini : minis) cout << mini.str() << ",";
    //            cout << endl;
    //        }

            for (Sequence mini : minis){


                auto tt0 = high_resolution_clock::now();
                //either the kmer has already been seen elsewhere or not
                if (index.find(mini) != index.end()){ //if it's not the first time this kmer is seen
                    auto ttt0 = high_resolution_clock::now();
                    long int place = index[mini];
                    re.new_mini(place, thread_id);

                    kmers[place].push_back(re.barcode); //each kmer contains the set of all barcodes containing this kmer

                    auto ttt1 = high_resolution_clock::now();
                    expansion_kmers_time += duration_cast<nanoseconds>(ttt1 - ttt0).count();

                }
                else {

                    auto ttt0 = high_resolution_clock::now();
                    mutex.lock();
                    index[mini] = kmers.size();
                    re.new_mini(kmers.size(), thread_id);
                    mutex.unlock();

                    kmers.push_back({re.barcode});
                    auto ttt1 = high_resolution_clock::now();
                    expansion_kmers_time += duration_cast<nanoseconds>(ttt1 - ttt0).count();
                }


                auto tt1 = high_resolution_clock::now();
                total_read_time += duration_cast<nanoseconds>(tt1 - tt0).count();

            }
            re.reset_minis_seq_thread(thread_id);


            count ++;
            if (count % (int(asize/100)+1) == 0 && verbose){
                cout << "Thread " << thread_id << " indexed " << double(count)/double(allreads.size())*100 << "% of its sequences" << endl;
            }
    }

    auto t1 = high_resolution_clock::now();

    float totk = 0;
    for (auto k : kmers){
        totk += std::set<long int>( k.begin(), k.end() ).size();
    }
    if(verbose){
        cout << "On average, each kmers is contained in " << totk/kmers.size() << " barcodes" << endl;
        cout << "In thread " << thread_id << " of index_reads, handling the memory of kmers took me reads " << total_read_time/1000000000 << "s out of " << duration_cast<microseconds>(t1 - t0).count()/1000000 << "s in total" <<  endl;
    }
}

//function eliminating all kmers that are abnormally often represented (more than twice the average). Also eliminate extra copies of barcodes present multiple times in one kmer
void convert_kmers(vector<vector<long int>>& res, vector<vector<long int>>& input, bool metagenomic, bool verbose){

    // a count to keep track of, on average, how many entries there are per kmers
    float coverage = 0;
    int count = 0;

    for (int l = 0 , ls = input.size() ; l < ls ; l++){

        set <long int> uniquek (input[l].begin(), input[l].end());
        res[l] = vector <long int> (uniquek.begin(), uniquek.end());

        int s = uniquek.size();
        if (s > 5){
            count++;
            coverage += uniquek.size();
        }
    }

    //now empty all the kmers that are anormally large (or else they will create false positives)
    if (!metagenomic){ //run this step only on non-metagenomic samples (on metagenomic samples it is to be expected that coverage may vary)
        float averageCoverage = coverage/count;
        int totrash = 0;
        for (int l = 0 , ls = res.size() ; l < ls ; l++){

            if (res[l].size() > 2*averageCoverage){
                res[l] = {};
                totrash ++;
            }

        }
        if(verbose){
            cout << "Eliminating for the classification all kmers present more than twice the average, i.e. " << 2*averageCoverage <<" times. Thread trashed " << totrash << " kmers out of " << count << " in total." << endl;
        }
    }
}
