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
void parse_reads(std::string fileReads, std::vector<std::vector<long long int>> &readClouds, std::vector <Read> &allreads, unordered_map <string, long int> &tagIDs, int min_length, int num_threads){

    //a count to inform the user how many reads will not be deconvolved (either because they do not have a tag or are too short)
    int totalNumberOfSequences = 0;
    auto t0 = high_resolution_clock::now();
    char format = '@'; //a character to keep track of whether the input file is a fasta or a fastq
    if ((fileReads.size()>6 && fileReads.substr(fileReads.size()-6,6) == ".fasta") || fileReads.substr(fileReads.size()-3,3) == ".fa"){
        format = '>';
    }

    ifstream in(fileReads);
    if (!in){
        cout << "problem reading files in index_reads, while trying to read " << fileReads << endl;
        throw std::invalid_argument( "Input file could not be read" );
     }

    long long int sequenceID = 0; //instead of storing the name of sequence in a string we're going to store them in long long int, it will be more efficient
    string lastnameofsequence; //used to see if it's the same as the current name of sequence: if it is, reads are paired

    string tag;

    long int tagID;

    string line;
    vector<string> buffer;

    while(getline(in, line)){

        if (line[0] == format){

            //then first we append the last read we saw

            if (to_deconvolve(buffer, format, min_length, tag)){ //to_deconvolve checks if the line should be deconvolved and if so returns the tag

                if (buffer[0] == lastnameofsequence){ //if it has the same name as the previous read, it means it is paired
                    allreads[allreads.size()-1].sequence_paired = buffer[1];
                    totalNumberOfSequences--;
                }
                else{

                    Read r(num_threads);

                    if (tagIDs.find(tag) == tagIDs.end()){ // if true, tag is not in the keys of tagIDs

                        tagIDs[tag] = readClouds.size(); //from now on this tag ID will be associated to this tag
                        tagID = readClouds.size();
                        readClouds.push_back({sequenceID});

                    }
                    else{
                        tagID = tagIDs[tag];
                        readClouds[tagID].push_back(sequenceID);
                    }

                    r.sequence = buffer[1];
                    r.barcode = tagID;
                    allreads.push_back(r);
                    sequenceID ++;

                    lastnameofsequence = buffer[0];

                }

            }

            //then we reset the buffer
            buffer = {line};

            totalNumberOfSequences ++;
        }
        else {
            buffer.push_back(line);
        }

    }

    //parse the last read
    if (to_deconvolve(buffer, format, min_length, tag)){

        if (buffer[0] == lastnameofsequence){ //if it has the same name as the previous read, it means it is paired
            allreads[allreads.size()-1].sequence_paired = buffer[1];
        }
        else{

            Read r(num_threads);

            if (tagIDs.find(tag) == tagIDs.end()){ // if true, tag is not in the keys of tagIDs

                tagIDs[tag] = readClouds.size(); //from now on this tag ID will be associated to this tag
                tagID = readClouds.size();
                readClouds.push_back({sequenceID});

            }
            else{
                tagID = tagIDs[tag];
                readClouds[tagID].push_back(sequenceID);
            }

            r.sequence = buffer[1];
            r.barcode = tagID;
            allreads.push_back(r);
            sequenceID ++;

            lastnameofsequence = buffer[0];
        }

    }

    auto t1 = high_resolution_clock::now();
    cout << "Finished parsing all the reads in " << duration_cast<seconds>(t1-t0).count() << "s, for a total of " << sequenceID << " reads to deconvolve. " << totalNumberOfSequences - sequenceID << " reads will not be deconvolved, either because they are too short or do not have a tag" << endl;
}

//function computing minimisers of some reads(one thread)
void compute_minimizers(int k, int h, int w, std::vector <Read> &allreads, int thread_id, int num_threads){

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

        if ((r-thread_id)%(100000*num_threads) == 0){
            cout << "Thread " << thread_id << " computed " << double(r)/double(allreads.size())*100 << "% of its minimizers" << endl;
        }
    }

    auto t1 = high_resolution_clock::now();

    cout << "Thread " << thread_id << " computed its minimizers in " << duration_cast<seconds>(t1-t0).count()<< "s" << endl;
}

//this function partially fills the index (kmers) with one thread
void index_kmers(int k, vector<vector<long int>> &kmers, vector <Read> &allreads, int thread_id, int num_threads){ //kmers is a vector of vector (and not vec of vec of vec) becuse we're focusing on 1 thread

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
            if (count % (int(asize/100)+1) == 0){
                cout << "Thread " << thread_id << " indexed " << double(count)/double(allreads.size())*100 << "% of its sequences" << endl;
            }
    }

    auto t1 = high_resolution_clock::now();

    float totk = 0;
    for (auto k : kmers){
        totk += std::set<long int>( k.begin(), k.end() ).size();
    }
    cout << "On average, each kmers is contained in " << totk/kmers.size() << " barcodes" << endl;

    cout << "In thread " << thread_id << " of index_reads, handling the memory of kmers took me reads " << total_read_time/1000000000 << "s out of " << duration_cast<microseconds>(t1 - t0).count()/1000000 << "s in total" <<  endl;
}

//function eliminating all kmers that are abnormally often represented (more than twice the average). Also eliminate extra copies of barcodes present multiple times in one kmer
void convert_kmers(vector<vector<long int>>& res, vector<vector<long int>>& input, bool metagenomic){

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
        cout << "Eliminating for the classification all kmers present more than twice the average, i.e. " << 2*averageCoverage <<" times. Thread trashed " << totrash << " kmers out of " << count << " in total." << endl;
    }
}


