#include "build_graph.h"
#include "cluster_graph.h"

#include <random> //for the iota function
#include <map>

using namespace std::chrono;

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::min;
using robin_hood::unordered_map;
using robin_hood::unordered_flat_map;
using std::map;
using std::string;
using std::array;
using std::unordered_set;
using std::set;

//function takes as input the reads and the index and deconvolves the barcodes the thread should deconvolve
void thread_deconvolve(short minCommonKmers, unordered_map <string, long int> &tagIDs, const vector <vector<long long int>> &readClouds, std::vector <Read> &reads, const vector<vector<vector<long int>>> &kmers, int thread_id, int num_thread, int dropout, string folderOut){

    int count = 0;
    for (auto p : tagIDs){

        //a condition because we want each thread to work separately
        if (p.second % num_thread == thread_id && readClouds[p.second].size() > dropout){
            vector <int> clusters (readClouds[p.second].size(), -1);
            build_graph(minCommonKmers, p.first, p.second, readClouds, reads, kmers, clusters, folderOut);
            //cout << "Deconvolved " << p.second << endl;

            if (count % 100 == 0) cout << "thread " << thread_id << " deconvolved " << count << " tags over " << readClouds.size() << " in total" << endl;
            count ++;
         }

    }
}

//the function takes as an input the list of all reads having the same tag, build the adjMatrix with the index and then cluters the graph: globally, it deconvolves one tag

void build_graph(short minCommonKmers, string tag, long int tagCloud, const vector <vector<long long int>> &readClouds, std::vector <Read> &reads, const vector<vector<vector<long int>>> &kmers, vector<int> &clusters, string folderOut){
	
	auto t0 = high_resolution_clock::now();

    int adjMatrixSize = clusters.size();
    vector<int> zeros (adjMatrixSize, 0);
    vector<vector<int>> adjMatrix (adjMatrixSize, zeros);
		
    unordered_map<long int, unordered_set<int>> matching_tags; //maps to each tagID the index of the reads of the cloud aligning with this tag
    build_adj_matrix(minCommonKmers, tagCloud, readClouds, reads, kmers, adjMatrix, matching_tags);

    auto t1 = high_resolution_clock::now();

    cluster_graph_chinese_whispers(adjMatrix, clusters, tag);
    //now store the result in the reads
    for (int r = 0 ; r < readClouds[tagCloud].size() ; r++ ){

        reads[readClouds[tagCloud][r]].barcode_extension = clusters[r];

    }

    auto t2 = high_resolution_clock::now();

 //   fast_clustering(tagCloud, readClouds, reads, kmers, clusters);

    auto t3 = high_resolution_clock::now();

//    cout << "Building adjacency matrix : " << duration_cast<microseconds>(t1-t0).count()/1000 << "ms, clustering the matrix : " << duration_cast<microseconds>(t2-t1).count()/1000 << "ms, fast clustering : " << duration_cast<microseconds>(t3-t2).count()/1000 << "ms" << endl;
	

    if (adjMatrix.size() > 3000){

        if (folderOut[folderOut.size()-1] != '/'){
            folderOut += '/';
        }

        string f = folderOut + "cluster_"+tag+"_adj.csv";
        string f2 = folderOut + "cluster_"+tag+"_nodes.csv";
        cout << "exporting..."  << adjMatrix.size()  << " "<< clusters.size()<< endl;

        export_as_CSV(adjMatrix, f, f2, clusters);
        f = folderOut + "cluster_"+tag+"_matching-tag.csv";
        export_as_CSV(matching_tags, f);
    }
}

//functions builds adjMatrix of one barcode, using the index
void build_adj_matrix(short minCommonKmers, long int tagCloud, const vector <vector<long long int>> &readClouds, std::vector <Read> &reads, const vector<vector<vector<long int>>> &kmers, vector<vector<int>> &adjMatrix, unordered_map<long int, unordered_set<int>>& matching_tags){

    auto t0 = high_resolution_clock::now();

    double d = 0;
    double d2 = 0;

    //int r = 0;
    for(int r = 0, sizer = readClouds[tagCloud].size(); r<sizer ; r++){


        unordered_map<long int, int> alreadySeen; //a map to keep track of how many times the read has already been attached to that tag: you need to have at least minCommonKmers common minimizer to attach a read to a tags
        long long int name = readClouds[tagCloud][r];

        vector<vector<long int>> &minis = reads[name].get_minis();

        for (short t = 0 ; t < kmers.size() ; t++){ //here, we iterate through all threads that built the index

            for (long int m : minis[t]){

                for (long int tag : kmers[t][m]){

                        auto tt0 = high_resolution_clock::now();
                        alreadySeen[tag] += 1;

                        if (alreadySeen[tag] == minCommonKmers){ // it would be equivalent to put >= here, but a bit slower
                            matching_tags[tag].emplace(r);
                        }
                        auto tt1 = high_resolution_clock::now();
                        d += duration_cast<nanoseconds>(tt1-tt0).count();
                }
            }
        }

    }
    matching_tags[tagCloud] = {}; //this line to avoid self-loops

    auto t1 = high_resolution_clock::now();

    //now build the interaction matrix between reads
    for (robin_hood::pair<const long int, unordered_set<int>> matchs : matching_tags){

//		cout << "\nLet's look at what matched with tag " << matchs.first << endl;
        int i = 0;
        for (int a  : matchs.second) {

            int j = 0;
            for (int b : matchs.second){

                if (j>i){
                    adjMatrix[a][b] += 1;
                    adjMatrix[b][a] += 1;
                }
                j++;
            }
            i++;
        }
    }
    auto t2 = high_resolution_clock::now();

    //cout << "While building adjMat, took me " << duration_cast<microseconds>(t1-t0).count() << "us to create matching tags, among it " << int(d/1000) << "us handling maps (against potentially "<< int(d2/1000) <<"us) " << duration_cast<microseconds>(t2-t1).count() << "us to build the adjMat " << endl;
}
