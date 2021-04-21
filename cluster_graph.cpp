#include "cluster_graph.h"
#include <random>
#include <set>

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::string;
using robin_hood::pair;
using std::set;
using robin_hood::unordered_map;

void cluster_graph(unordered_map<long int, std::set<int>> &matching_tags, vector<int> &clusters){
	
	int adjMatrixSize = clusters.size();
	vector<int> zeros (adjMatrixSize, 0);
	vector<vector<int>> adjMatrix (adjMatrixSize, zeros);
	
    int reject_tag_threshold = 30000000; //if a readCloud overlaps with too many reads in your readCloud, consider it a fluke
	vector<int> strengths_of_links = {0}; //strengths_of_links[p] contains how many links are at least as strong as p
	
	
	//building the interaction matrix between reads
    for (pair<const long int, set<int>> matchs : matching_tags){
		
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
    int link_between_reads_threshold = std::max(1,find_threshold(0.8, strengths_of_links)); //if two reads share at least this number of common overlapping tags, consider them linked
//	cout << "Threshold : " << link_between_reads_threshold << endl;
	
    //int link_between_reads_threshold = 1;
	
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
	
    if (adjMatrix.size()>20){
        string id = std::to_string(int(rand()%30));
        string f = "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/cluster_"+id+"_adj.csv";
        export_as_CSV(adjMatrix, f);
        f = "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/cluster_"+id+"_matching-tag.csv";
        export_as_CSV(matching_tags, f);
    }
	//export_as_CSV(adjMatrix, "evalResultGraphs/cluster.csv");
		
	find_connected_components(adjMatrix, clusters);
	

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

void cluster_graph_chinese_whispers(unordered_map<long int, std::set<int>> &matching_tags, vector<int> &clusters){

    int adjMatrixSize = clusters.size();
    vector<int> zeros (adjMatrixSize, 0);
    vector<vector<int>> adjMatrix (adjMatrixSize, zeros);

    //building the interaction matrix between reads
    for (pair<const long int, set<int>> matchs : matching_tags){

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

    // this can be suppressed, it is just for exporting examples
    for (int i =0 ; i<adjMatrix.size(); i++){
        for (int j = 0 ; j<adjMatrix[i].size() ; j++){
            if (adjMatrix[i][j] < 2){
                adjMatrix[i][j] = 0;
            }
        }
    }


    float nb_of_nodes = clusters.size();
    float nb_of_nodes_changing = nb_of_nodes; //an index keeping count of how many nodes changed in the last iteration, to estimate when convergence is reached
    int nb_of_iterations = 0; //another index keeping count of how many iterations there was, to be sure not to get trapped in a oscillating result

    vector<int> order (clusters.size()); //a vector telling at each step in which order to go through the nodes
    std::iota(order.begin(), order.end(), 0); //fill this vector with increasing int starting from 0 (we'll shuffle that later)

    std::iota(clusters.begin(), clusters.end(), 0); //each node begins as its own class

    while (nb_of_iterations < 100 && nb_of_nodes_changing > 0.05*nb_of_nodes){

        nb_of_iterations++;
        nb_of_nodes_changing = 0; //it will then be increased through this loop

        std::random_shuffle(order.begin(), order.end()); //define the random order with which to look at the node

        for (int node : order){ //go through the nodes in a random order

            //compute the frequencies of the clusters in the vicinity
            vector <int> frequency_neighbors(nb_of_nodes, 0);

            for (int neighbor = 0 ; neighbor < nb_of_nodes ; neighbor++){
                //if (adjMatrix[node][neighbor] > 1){
                frequency_neighbors[clusters[neighbor]] += adjMatrix[node][neighbor];
                //}

            }

            int occurences_of_most_frequent_neighbor = 2; //do not start at 0 is equivalent to saying that two sequences share a link if and only if they share more than 0 kmers
            vector <int> most_frequent_neighbors;

            //find the most frequent clusters in the vicinity by going through frequency_neighbor
            for (int neighbor = 0 ; neighbor < nb_of_nodes ; neighbor++){

                if (frequency_neighbors[neighbor] > occurences_of_most_frequent_neighbor){
                    occurences_of_most_frequent_neighbor = frequency_neighbors[neighbor];
                    most_frequent_neighbors = {neighbor};
                }
                else if (frequency_neighbors[neighbor] == occurences_of_most_frequent_neighbor){
                    most_frequent_neighbors.push_back(neighbor);
                }

            }
            //randomly pick one of the most frequent neighbor
            if (most_frequent_neighbors.size() > 0){
                int newCluster = most_frequent_neighbors[rand()%most_frequent_neighbors.size()];
                //cout << "Next cluster : " << newCluster << ", old one : " << clusters[node] << endl;
                if (newCluster != clusters[node]){
                    clusters[node] = newCluster;
                    nb_of_nodes_changing++;
                }
            }
        }

    }


    //cout << "Exporting..." << endl;
//    if (adjMatrix.size()>20){
//        string id = std::to_string(int(rand()%30));
//        string f = "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/cluster_"+id+"_adj.csv";
//        export_as_CSV(adjMatrix, f);
//        f = "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/cluster_"+id+"_matching-tag.csv";
//        export_as_CSV(matching_tags, f);
//    }

}
