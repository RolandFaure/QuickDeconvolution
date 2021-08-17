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
    //measure_graph_building_time(20,3,40, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/eval/reads_1Mb_cov25_redundance4.fastq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/",  "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );
    //measure_graph_building_time(20,3,40, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/datasets/H_numata/barcoded.tiny.fastq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );
    //measure_graph_building_time(20,3,40, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/datasets/mock_metagenomes/10M.data1_atgctgaaq.small.fq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );

    int k = 20 , w = 40 , h = 3, t = 1, a=0;
    string infile, outfolder, outfile;
    auto cli = (
            required("-i", "--input-file") & opt_value("i", infile),
            required("-f", "--output-folder") & opt_value("f", outfolder),
            required("-o", "--output-file") & opt_value("o", outfile),
            option("-k", "--kmers-length").doc("size of kmers") & opt_value("k", k),
            option("-w", "--window-size").doc("size of window guaranteed to contain at least one minimizing kmer") & opt_value("w", w),
            option("-d", "--density").doc("on average 1/2^d kmers are sparse kmers") & opt_value("d", h),
            option("-t", "--threads").doc("number of threads") & opt_value("t", t),
            option("-a", "--dropout").doc("QD does not try to deconvolve clouds smaller than this value [default:0]") & opt_value("a", a)
        );

    if(!parse(argc, argv, cli)) {
        cout << "Could not parse the arguments" << endl;
        cout << make_man_page(cli, argv[0]);
    }
    else {
        cout << "Launching deconvolution, with arguments : k=" <<k << " d=" << h << " w=" << w << " t=" << t << " infile=" << infile << " outfolder=" << outfolder << " outfile=" << outfile << endl;
        measure_graph_building_time(k, h, w, t,a, infile, outfolder, outfile);
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
