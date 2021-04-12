#include "cluster_graph.h"

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::string;
using robin_hood::pair;
using robin_hood::unordered_map;

void cluster_graph(unordered_map<long int, std::list<int>> &matching_tags, vector<int> &clusters){
	
	int adjMatrixSize = clusters.size();
	vector<int> zeros (adjMatrixSize, 0);
	vector<vector<int>> adjMatrix (adjMatrixSize, zeros);
	
	int reject_tag_threshold = 30; //if a readCloud overlaps with too many reads in your readCloud, consider it a fluke
	vector<int> strengths_of_links = {0}; //strengths_of_links[p] contains how many links are at least as strong as p
	
	
	//building the interaction matrix between reads
	for (pair<long int, list<int>> matchs : matching_tags){
		
		int mss = matchs.second.size();
		if (mss < reject_tag_threshold){
	//		cout << "\nLet's look at what matched with tag " << matchs.first << endl;
			int i = 0;
			for (int a  : matchs.second) {
				
				int j = 0;
				for (int b : matchs.second){
					
					if (j>i){
						adjMatrix[a][b] += 1;
						adjMatrix[b][a] += 1;
						if (adjMatrix[a][b]<strengths_of_links.size()){
							strengths_of_links[adjMatrix[a][b]] += 1;
						}
						else{
							strengths_of_links.push_back(1);
						}
					}
					j++;
				}
				i++;
			}
		}
	}
	
//	for (int i : strengths_of_links) cout << i << " ";
	int link_between_reads_threshold = find_threshold(0.5, strengths_of_links); //if two reads share at least this number of common overlapping tags, consider them linked
//	cout << "Threshold : " << link_between_reads_threshold << endl;
	
	//link_between_reads_threshold = 3;
	
	//building the adjacency matrix : if enough common tags, the reads are considered linked
	//link together the reads with enough common barcodes
	for (int i = 0 ; i < adjMatrix.size() ; i++){
				
		for (int j = 0 ; j < adjMatrix.size() ; j++){
			if (adjMatrix[i][j] >= link_between_reads_threshold){
				adjMatrix[i][j] = 1;
			}
			else{
				adjMatrix[i][j] = 0;
			}
		}
	}
	
	//export_as_CSV(adjMatrix, "evalResultGraphs/cluster.csv");
		
	find_connected_components(adjMatrix, clusters);
	
//	string f = "evalResultGraphs/"+to_string(t)+"_adj.csv";
//	export_as_CSV(adjMatrix, f);
}

//a function to find connected components in an adjacency matrix using b
void find_connected_components(vector<vector<int>> &adjMatrix, vector<int> &clusters){
	
	int ncomponent = 0;
	for (int i = 0, sizei = adjMatrix.size() ; i < sizei ; i++){
		if (clusters[i] == -1){
			clusters[i] = ncomponent;
			extend_cluster_from_seed(i, adjMatrix, clusters);
			ncomponent += 1;
		}
	}
	
}

//recursively extend a given connected component from the node seed
void extend_cluster_from_seed(int seed, vector<vector<int>> &adjMatrix, vector<int> &clusters){
	
	for (int neighbor, sizen = adjMatrix.size() ; neighbor < sizen ; neighbor ++){
		if (adjMatrix[seed][neighbor] == 1 and clusters[neighbor] != clusters[seed]){
			
//			if (clusters[neighbor] != -1)cout << "Problem in recursive component finding " << endl;
			
			clusters[neighbor] = clusters[seed];
			extend_cluster_from_seed(neighbor, adjMatrix, clusters);
		}
	}
	
}

//gives what threshold should we use to keep the given proportion of best links 
int find_threshold(float proportion, vector<int> &strengths_of_links){
	
	int total_number = strengths_of_links[2]; //all links are at least 2 strong, since links 1 strong are not not trustworthy
	int kept_number = int(total_number*proportion);
	
	//proceed by dichotomy to find the wanted threshold
	int a = 1;
	int b = strengths_of_links.size()-1;
	
	while (b-a > 1){
		int c = int((a+b)/2);
		if (strengths_of_links[c]>kept_number){
			a = c;
		}
		else{
			b = c;
		}
	}
	
	return b;
}
