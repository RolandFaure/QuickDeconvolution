#include "tests.h"

using std::cout;
using std::endl;
using std::vector;
using std::list;
using robin_hood::pair;
using std::string;
using std::array;
using robin_hood::unordered_map;
using std::set;

vector<vector<int>> true_adjMatrix(vector<long long int> cloud, vector <Read> allReads){
	
//	vector<int> zeros (cloud.size(), 0);
//	vector<vector<int>> adjMatrix (cloud.size(), zeros);
	
//	unordered_map<string, vector<int>> matching_tags;
	
//	int index = 0;
//	for (long long int name : cloud){
//        long int t = allReads[name].barcode;
//		matching_tags[t].push_back(index);
//		index ++;
//	}
	
//    for (pair<const string, vector<int>> matchs : matching_tags){
		
//		//cout << "Let's look at what matched with tag " << matchs.first << endl;
//		for (int i = 0; i<matchs.second.size()-1 ; i++) {
//			int a = matchs.second[i];

//			for (int j = i+1 ; j < matchs.second.size() ; j++){
//				int b = matchs.second[j];
//				//cout << reads[readCloud[matchs.second[j]]] << endl;
//				adjMatrix[a][b] = 1;
//				adjMatrix[b][a] = 1;
//			}
//		}
//	}
//	return adjMatrix;
}

void rapid_check(){
	
//	char alphabet[4] = {'A', 'C', 'G', 'T'};
//	chromosomes_to_file(100000, alphabet, 4, 1, "eval/genome_smallTest.fasta");
//	draw_fragments(150, 400, 5, 0.2, 15000, 20000, 2, "eval/genome_smallTest.fasta", "eval/reads_smallTest");
	
    vector<set<long int>> kmers;
	vector <vector<long long int>> readClouds;
	vector <Read> allreads;
	int k = 20;
	int w = 10;
    int h = 3;
	
    index_reads(k, h, w, "eval/reads_10Mb_cov25_redundance4.fasta", kmers, readClouds, allreads);
	cout << "Finished indexing" << endl;
	
	long int index = 0;
	for (vector<long long int> cloud : readClouds){
		
		if (index == 0 /*&& cloud.size() < 120*/ ){
            vector<int> clusters(cloud.size(), -1);
            build_graph(3, index, cloud, allreads, kmers, clusters);
			
			vector<vector<int>> adjMatrix_t = true_adjMatrix(cloud, allreads);
		
			cout << "Computed adj_matrix of tag " << index << " : " << endl;
//			for (int i=0 ; i < adjMatrix.size() ; i++){
//				for (int j = 0 ; j < adjMatrix.size() ; j++){
//					cout << adjMatrix[i][j] << "\t";
//				}
//				cout << endl;
//			}
//			export_as_SIF(adjMatrix, "evalResultGraphs/"+to_string(index)+".sif");
			
//			cout << "True adjacent matrix of the tag : " << endl;
//			for (int i=0 ; i < adjMatrix_t.size() ; i++){
//				for (int j = 0 ; j < adjMatrix_t.size() ; j++){
//					cout << adjMatrix_t[i][j] << "\t";
//				}
//				cout << endl;
//			}
            export_as_CSV(adjMatrix_t, "evalResultGraphs/"+std::to_string(index)+"_true.csv");
			cout<<"Done for tag " << index << endl;
		}
		index ++;

	}
}

