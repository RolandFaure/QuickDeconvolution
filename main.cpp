#include <stdio.h>
#include "create_experiments.h"
#include "measure_speed.h"
#include "generate_random_sequence.h"
#include "draw_10x.h"
#include "tests.h"
#include "cluster_graph.h"
#include "clipp.h" //library to build command line interfaces

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::pair;
using std::string;
using robin_hood::unordered_map;
using namespace clipp;

int main(int argc, char *argv[])
{

    //create_exps();
    //rapid_check();

    int num_threads = 2;
    //measure_graph_building_time(20,4,40,50, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/eval/reads_10Mb_cov50_redundance4.fastq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/",  "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );
    measure_graph_building_time(20,3,40,50, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/datasets/mock_metagenomes/mock6_lsq.R1.short.fastq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );
    //measure_graph_building_time(20,3,40,50, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/datasets/trick_reads.fastq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );

    int k = 20 , w = 40 , h = 4, c = 50, t = 1;
    string infile, outfolder, outfile;
    auto cli = (
            value("input file to deconvolve", infile),
            value("deconvolved output folder", outfolder),
            value("deconvolved output file", outfile),
            option("-k").set(k).doc("size of kmers"),
            option("-w").set(w).doc("size of window guaranteed to contain at least one minimizing kmer"),
            option("-h").set(h).doc("on average 1/2^ĥ kmers are minimizing kmers"),
            option("-t").set(t).doc("number of threads")
        );

    if(!parse(argc, argv, cli)) cout << make_man_page(cli, argv[0]);
    else {
        measure_graph_building_time(k, h, w, c, t, infile, outfolder, outfile);
    }


//    if (argc > 4){
//        cout << "Launching deconvolution, with arguments : 20, 4, 40, 50," << std::stoi(argv[4]) << ", " << argv[1] << ", " << argv[2] << ", " << argv[3] << endl;
//        measure_graph_building_time(20, 4, 40, 50, std::stoi(argv[4]), argv[1], argv[2], argv[3]);
//    }
//    else cout << "You must give as an argument the fastq file I will try to deconvolve, the path to output folder, path to output file and the number of threads" << endl;
    //systematic_times(30);
    //draw_fragments(150, 400, 25, 0.2, 20000, 5000, 4, "eval/genome_10Mb.fasta", "eval/reads_10Mb_cov25_redundance4.fasta");
    return 0;
}
