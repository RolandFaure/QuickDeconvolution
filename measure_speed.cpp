#include "measure_speed.h"

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::string;
using std::array;
using std::set;
using robin_hood::pair;
using robin_hood::unordered_map;
using std::thread;
using std::ref; //to pass references in threads
using std::this_thread::sleep_for; //to pause the program
using namespace std::chrono;



float measure_graph_building_time(int k, int h, int w, int c, int num_threads, string readsFile, string folderOut){

    double timeGraph = 0;
    auto t0 = high_resolution_clock::now();

    vector<vector<long int>> kmers;
    vector <vector<long long int>> readClouds;
    vector <Read> allreads;

    unordered_map <string, long int> tagIDs;
    index_reads(k, h, w, readsFile, kmers, readClouds, allreads, tagIDs);


    //for now kmers is a vector of sets, but iterating through sets is not very fast, so we'll convert that into a vector of vectors
    cout << "Converting kmers to make it faster to iterate" << endl;
    vector<vector<long int>> kmersV(kmers.size());
    double meansize = 0;
    for (int l = 0 , ls = kmers.size() ; l < ls ; l++){

        set <long int> uniquek (kmers[l].begin(), kmers[l].end());
        vector <long int> kv (uniquek.begin(), uniquek.end());

        if (kv.size() > c){
            kv.erase(kv.begin()+c , kv.end()); //only keep the first c elements
        }
        meansize += kv.size();
        kmersV[l] = kv;
    }

    vector<vector<long int>> ().swap(kmers); //this frees up the memory taken by kmers, since we'll only be using kmersV from now on
    cout << "Finished converting kmers, mean size of kmers : " << meansize/kmersV.size() << endl;


	auto t1 = high_resolution_clock::now();

    vector<thread> threads;

    for (int i = 0 ; i < num_threads ; i++){

        threads.push_back(thread(thread_deconvolve, 3, ref(tagIDs), ref(readClouds), ref(allreads), ref(kmersV), i, num_threads, folderOut));

    }

    //now join all the threads
    for (vector<thread>::iterator it = threads.begin() ; it != threads.end() ; ++it)
    {
        it->join();
    }
	
	auto t2 = high_resolution_clock::now();
	
    cout << "Indexation time : " << duration_cast<seconds>(t1 - t0).count() << "s, alignment time " << duration_cast<milliseconds>(t2 - t1).count()/1000 << "s " << endl;
	return duration_cast<seconds>(t2 - t0).count();

}

void systematic_times(int k, int w){

    string exp1 = "eval/reads_1Mb_cov25_redundance4_1.fasta";
    string exp2 = "eval/reads_1Mb_cov50_redundance4_1.fasta";
    string exp3 = "eval/reads_1Mb_cov50_redundance15_1.fasta";

    string exp4 = "eval/reads_10Mb_cov25_redundance4_1.fasta";
    string exp5 = "eval/reads_10Mb_cov50_redundance4_1.fasta";
    string exp6 = "eval/reads_10Mb_cov50_redundance15_1.fasta";

    string exp7 = "eval/reads_100Mb_cov25_redundance4_1.fasta";
    string exp8 = "eval/reads_100Mb_cov50_redundance4_1.fasta";
    string exp9 = "eval/reads_100Mb_cov50_redundance15_1.fasta";

    string exp10 = "eval/reads_1Gb_cov25_redundance4_1.fasta";
    string exp11 = "eval/reads_1Gb_cov50_redundance4_1.fasta";
    string exp12 = "eval/reads_1Gb_cov50_redundance15_1.fasta";

//	cout << exp1 << "\t" << measure_graph_building_time(k, w, exp1) << "s" << endl;
//	cout << exp2 << "\t" << measure_graph_building_time(k, w, exp2) << "s" << endl;
//	cout << exp3 << "\t" << measure_graph_building_time(k, w, exp3) << "s" << endl;

//	cout << exp4 << "\t" << measure_graph_building_time(k, w, exp4) << "s" << endl;
//	cout << exp5 << "\t" << measure_graph_building_time(k, w, exp5) << "s" << endl;
//	cout << exp6 << "\t" << measure_graph_building_time(k, w, exp6) << "s" << endl;

//	cout << exp7 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	cout << exp8 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	cout << exp9 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//
//	cout << exp10 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	cout << exp11 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	cout << exp12 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
}
