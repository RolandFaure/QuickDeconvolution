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
//using std::atomic; //to use atomic types (multithreading)



float measure_graph_building_time(int k, int h, int w, int c, int num_threads, string readsFile, string folderOut, string fileOut){

    double timeGraph = 0;
    auto t0 = high_resolution_clock::now();

    vector<vector<vector<long int>>> kmers (num_threads); //different index, one for each thread
    vector <vector<long long int>> readClouds;
    vector <Read> allreads;

    unordered_map <string, long int> tagIDs;

    //first parse all the reads from the file
    int min_length = w+k; //then minimum length of a read to have it deconvolved
    parse_reads(readsFile, readClouds, allreads, tagIDs, min_length, num_threads);


    //then compute all the minimizers (parallely)

    auto t0_5 = high_resolution_clock::now();
    vector<thread> threadsMini;

    for (int i = 1 ; i < num_threads ; i++){

        threadsMini.push_back(thread(compute_minimizers, k, h, w, ref(allreads), i, num_threads));

    }
    compute_minimizers(k,h,w, ref(allreads), 0, num_threads); //also use the main thread, no reason to wait

        //now join all the threads
    for (vector<thread>::iterator it = threadsMini.begin() ; it != threadsMini.end() ; ++it)
    {
        it->join();
    }

    auto t1 = high_resolution_clock::now();

    //finish by building the index (kmers) (parallely)

    vector<thread> threadsIdx;

    for (int i = 1 ; i < num_threads ; i++){

        threadsIdx.push_back(thread(index_kmers,k, ref(kmers[i]), ref(allreads), i, num_threads));

    }
    index_kmers(k, kmers[0], allreads, 0, num_threads);


        //now join all the threads
    for (vector<thread>::iterator it = threadsIdx.begin() ; it != threadsIdx.end() ; ++it)
    {
        it->join();
    }

    cout << "First, allreads minis : " << allreads[0].get_minis()[0].size()<< endl;

    auto t1_5 = high_resolution_clock::now();

    cout << "In total, indexing reads took me " << duration_cast<seconds>(t1_5-t0).count() << "s, among which " << duration_cast<seconds>(t0_5-t0).count() << "s for parsing, " << duration_cast<seconds>(t1-t0_5).count() << "s for finding minimizers and " << duration_cast<seconds>(t1_5-t1).count() << "s for putting all that in an index"  << endl;

    //for now kmers[i][j] may contain replicate (in repetitive regions), get rid of them to iterate much faster

    vector<thread> threadsConvert;
    vector<vector<vector<long int>>> kmersV(kmers.size());

    for (int t = 1 ; t < num_threads ; t++){

        kmersV[t] = vector<vector<long int>> (kmers[t].size());
        threadsConvert.push_back(thread(convert_kmers,ref(kmersV[t]), ref(kmers[t]) , c));

    }
    kmersV[0] = vector<vector<long int>> (kmers[0].size());
    convert_kmers(kmersV[0], kmers[0], c);

        //now join all the threads
    for (vector<thread>::iterator it = threadsConvert.begin() ; it != threadsConvert.end() ; ++it)
    {
        it->join();
    }


    vector<vector<vector<long int>>> ().swap(kmers); //this frees up the memory taken by kmers, since we'll only be using kmersV from now on


    auto t10 = high_resolution_clock::now();

    vector<thread> threads;

    for (int i = 1 ; i < num_threads ; i++){

        threads.push_back(thread(thread_deconvolve, 3, ref(tagIDs), ref(readClouds), ref(allreads), ref(kmersV), i, num_threads, folderOut));

    }
    thread_deconvolve(3, ref(tagIDs), ref(readClouds), ref(allreads), ref(kmersV), 0, num_threads, folderOut);

    //now join all the threads
    for (vector<thread>::iterator it = threads.begin() ; it != threads.end() ; ++it)
    {
        it->join();
    }
	
	auto t2 = high_resolution_clock::now();

    //now write the output
    cout << "Finished deconvolving, now writing to the output" << endl;
    output(readsFile, fileOut, allreads);

    auto t3 = high_resolution_clock::now();
	
    cout << "Indexation time : " << duration_cast<seconds>(t1_5 - t0).count() << "s, conversion time "<< duration_cast<seconds>(t10 - t1_5).count() <<"s, alignment time " << duration_cast<milliseconds>(t2 - t10).count()/1000 << "s, output time " << duration_cast<milliseconds>(t3 - t2).count()/1000 << "s" << endl;
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
