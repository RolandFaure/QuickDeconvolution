#include "measure_speed.h"

using namespace std;
using namespace std::chrono; 

float measure_graph_building_time(int k, int w, string readsFile){
	
	double timeGraph = 0;
	auto t0 = high_resolution_clock::now();
	
	unordered_map <Sequence, vector<vector<Hit>>, Sequence::HashFunction> idx;
	vector <vector<long long int>> readClouds;
	vector <Read> allreads;
		
	index_reads(k, w, readsFile, idx, readClouds, allreads);
	
	auto t1 = high_resolution_clock::now();

	long int index = 0;
	for (vector<long long int> cloud : readClouds){
		
		//if (index == 0){
			
		auto tt1 = high_resolution_clock::now();
		vector<int> clusters = build_graph(k, w, index, cloud, allreads, idx);
		auto tt2 = high_resolution_clock::now();
		timeGraph += duration_cast<microseconds>(tt2 - tt1).count();
//		cout << "Treating tag number " << index << endl;
//			for (int i=0 ; i<adjMatrix.size() ; i++){for(int j=0 ; j<adjMatrix.size() ; j++){cout << adjMatrix[i][j] << "\t";} cout << endl;}
//		cout << endl;
		
		index ++;
		if (index%100 == 0){
			cout << "Treated " << index << " tags" << endl;
		}
		//}
	}
	
	auto t2 = high_resolution_clock::now();
	
	cout << "Indexation time : " << duration_cast<seconds>(t1 - t0).count() << "s, alignment time " << duration_cast<seconds>(t2 - t1).count() << "s over a total of " << index << " tags, taking on average " << timeGraph/index << "us to build one graph" << endl;
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
	
	cout << exp1 << "\t" << measure_graph_building_time(k, w, exp1) << "s" << endl;
	cout << exp2 << "\t" << measure_graph_building_time(k, w, exp2) << "s" << endl;
	cout << exp3 << "\t" << measure_graph_building_time(k, w, exp3) << "s" << endl;
	
	cout << exp4 << "\t" << measure_graph_building_time(k, w, exp4) << "s" << endl;
	cout << exp5 << "\t" << measure_graph_building_time(k, w, exp5) << "s" << endl;
	cout << exp6 << "\t" << measure_graph_building_time(k, w, exp6) << "s" << endl;
	
//	cout << exp7 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	cout << exp8 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	cout << exp9 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	
//	cout << exp10 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	cout << exp11 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
//	cout << exp12 << "\t" << measure_graph_building_time(k, exp1) << "s" << endl;
}
