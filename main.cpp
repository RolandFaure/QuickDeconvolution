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

    //create_exps();
    //rapid_check();

    //measure_graph_building_time(20,4,40,50, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/eval/reads_1Mb_cov25_redundance4.fastq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/");
    if (argc > 2){
        measure_graph_building_time(20, 4, 40, 50, argv[1], argv[2]);

    }
    else cout << "You must give as an argument the fastq file I will try to deconvolve and the path to output files" << endl;
    //systematic_times(30);
    //draw_fragments(150, 400, 25, 0.2, 20000, 5000, 4, "eval/genome_10Mb.fasta", "eval/reads_10Mb_cov25_redundance4.fasta");
    return 0;
}
