#include "build_graph.h"
#include "cluster_graph.h"

//#include <lemon/smart_graph.h> //graph library
//#include <lemon/list_graph.h>

using namespace std::chrono;

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::min;
using robin_hood::unordered_map;
using std::string;
using std::array;
using std::set;
//using namespace lemon;

//the function takes as an input the list of all reads having the same tag
vector<int> build_graph(int k, int w, long int tagCloud, const std::vector<long long int>& readCloud, const std::vector <Read> &reads, vector<vector<long int>> &kmers, vector<int> &clusters){
	
	long int mini_time = 0;
	auto t0 = high_resolution_clock::now();
		
    unordered_map<long int, set<int>> matching_tags; //maps to each tagID the index of the reads of the cloud aligning with this tag
	
	float total_read_time = 0;
	
	//int r = 0;
	for(int r = 0, sizer = readCloud.size(); r<sizer ; r++){

        unordered_map<long int, int> alreadySeen; //a map to keep track of how many times the read has already been attached to that tag: you need to have at least 3 common minimizer for estimating there is an overlap
		long long int name = readCloud[r];
		
        for (long int m : reads[name].minis){

            for (long int tag : kmers[m]){

                alreadySeen[tag] += 1;

                if (alreadySeen[tag] > 2){
                    matching_tags[tag].emplace(r);
                }
            }

        }

    }
    matching_tags[tagCloud] = {}; //this line to avoid self-loops
	
    cluster_graph_chinese_whispers(matching_tags, clusters);
	
//	int n = 0;
////	cout << "sequence of read 0 : " << fullnum2str(reads[0].sequence) << endl;
////	cout << "sequence of read 0' : " << fullnum2str(reads[1].sequence) << endl;
//	if (tagCloud == 0){
//		for (pair<long int, list<int>> matchs : matching_tags){
//			for (int i: matchs.second){
//				if (readCloud[i] == 0){
//					cout << "Read 0 matched with tag " << matchs.first << endl;
//					n++;
//				}
//				if (readCloud[i] == 1){
//					cout << "Read 0' matched with tag " << matchs.first << endl;
//					n++;
//				}
//			}
//		}
//	}
//	cout << "reads 0 matched with " << n << " tags" << endl;

	auto t2 = high_resolution_clock::now();
	
	//cout << "push_back time : " << total_read_time << endl;
    //cout << "Alignement time : " << duration_cast<microseconds>(t1 - t0).count() <<"us, total read accession : " << total_read_time << ", total minimizer time : " << mini_time << endl;
	return clusters;
}
