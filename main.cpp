#include <stdio.h>
#include "create_experiments.h"
#include "measure_speed.h"
#include "generate_random_sequence.h"
#include "draw_10x.h"
#include "tests.h"
#include "cluster_graph.h"

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::pair;
using std::string;
using robin_hood::unordered_map;

int main(int argc, char *argv[])
{
//	vector<vector<int>> adjMatrix = {{0,1,0,0},{1,0,0,1},{0,0,0,1},{0,1,1,0}};
//	vector<int> clusters = {-1,-1,-1,-1};
//	find_connected_components(adjMatrix, clusters);
//	for (int i = 0 ; i < clusters.size() ; i++) {
//		cout << clusters[i] << "\t";
//	}
//	cout << endl;

    //create_exps();
    //rapid_check();

    measure_graph_building_time(20,4,40,50, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/eval/reads_10Mb_cov25_redundance4.fasta");
    if (argc > 1){
        //measure_graph_building_time(20, 4, 40, argv[1]);

    }
    else cout << "You must give as an argument the fastq file I will try to deconvolve" << endl;
    //systematic_times(30);
    //draw_fragments(150, 400, 25, 0.2, 20000, 5000, 4, "eval/genome_10Mb.fasta", "eval/reads_10Mb_cov25_redundance4.fasta");
    return 0;
}
