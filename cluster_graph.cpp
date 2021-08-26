#include "cluster_graph.h"
#include <random>
#include <set>

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::string;
using robin_hood::pair;
using std::unordered_set;
using robin_hood::unordered_map;

using namespace std::chrono;

//takes an adjacency matrix and clusters it with the chinese-whispers algorithm
void cluster_graph_chinese_whispers(vector<vector<int>> &adjMatrix, vector<int> &clusters, string &tag){

    long long unsigned int bas = 0;

    //now we're done, let's whisper
    float nb_of_nodes = clusters.size();
    float nb_of_nodes_changing = nb_of_nodes; //an index keeping count of how many nodes changed in the last iteration, to estimate when convergence is reached
    int nb_of_iterations = 0; //another index keeping count of how many iterations there was, to be sure not to get trapped in a oscillating result

    vector<int> order (clusters.size()); //a vector telling at each step in which order to go through the nodes
    std::iota(order.begin(), order.end(), 0); //fill this vector with increasing int starting from 0 (we'll shuffle that later)

    std::iota(clusters.begin(), clusters.end(), 0); //each node begins as its own class

    while (nb_of_iterations < 100 && nb_of_nodes_changing > 0.05*nb_of_nodes){

        nb_of_iterations++;
        nb_of_nodes_changing = 0; //it will then be increased through this loop

        auto t0 = high_resolution_clock::now();

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
        auto t1 = high_resolution_clock::now();

        bas += duration_cast<microseconds>(t1-t0).count();

    }

    //now the clusters have very high id (like 639, 13 and 1008) -> let's renumber them 1, 2 and 3
    //also count the size of the clusters and consider too small clusters as non-identified reads and attach the barcode 0
    unordered_map<int, pair<int,int>> indices;
    short clustID = 1; //start with 1 because ID 0 is dedicated to reads we do not know how to cluster
    for (int i = 0 ; i < clusters.size() ; i++){
        if (indices.find(clusters[i]) == indices.end()){
            indices[clusters[i]] = {clustID,1}; //the second memeber of the pair represents the number of reads in the cluster
            clustID++;
        }
        clusters[i] = indices[clusters[i]].first;
        indices[clusters[i]].second += 1;
    }
    for (int i = 0 ; i < clusters.size() ; i++){
        if (indices[clusters[i]].second <= 2){
            clusters[i] = 0;
        }
    }

    //cout << "clustered in " << nb_of_iterations << " iterations, basic time : " << bas/1000 << " ms" << endl;

}
