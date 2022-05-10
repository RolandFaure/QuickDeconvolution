#include <stdio.h>
#include "measure_speed.h"
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

    //int num_threads = 2;
    //measure_graph_building_time(20,3,40, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/eval/reads_1Mb_cov25_redundance4.fastq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/",  "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );
    //measure_graph_building_time(20,3,40, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/datasets/H_numata/barcoded.tiny.fastq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );
    //measure_graph_building_time(20,3,40, num_threads,"/home/zaltabar/Documents/Ecole/X/4A/stage_M2/datasets/mock_metagenomes/10M.data1_atgctgaaq.small.fq", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/evalGraphs/output.tsv" );

    int k = 20 , w = 40 , h = 3, t = 1, a=0;
    bool metagenomic = false;
    string infile, outfolder, outfile;
    auto cli = (
            required("-i", "--input-file").doc("input file (mandatory)") & opt_value("i", infile),
            required("-o", "--output-file").doc("file to write the output (mandatory)") & opt_value("o", outfile),
            option("-k", "--kmers-length").doc("size of kmers [default:20]") & opt_value("k", k),
            option("-w", "--window-size").doc("size of window guaranteed to contain at least one minimizing kmer [default:40]") & opt_value("w", w),
            option("-d", "--density").doc("on average 1/2^d kmers are sparse kmers [default:3]") & opt_value("d", h),
            option("-t", "--threads").doc("number of threads [default:1]") & opt_value("t", t),
            option("-a", "--dropout").doc("QD does not try to deconvolve clouds smaller than this value [default:0]") & opt_value("a", a),
            option("-m", "--metagenome").set(metagenomic).doc("Use this option on metagenomic samples")
        );

    if(!parse(argc, argv, cli)) {
        cout << "Could not parse the arguments" << endl;
        cout << make_man_page(cli, argv[0]);
    }
    else {
        //now deduce the output folder from the output file
        bool write = false;
        for (int i = 1 ; i <= outfile.size() ; i++){
            if (outfile[outfile.size()-i] == '/'){
                write = true;
            }
            if (write){
                outfolder = outfile[outfile.size()-i] + outfolder;
            }
        }

        cout << "Launching deconvolution, with arguments : k=" <<k << " d=" << h << " w=" << w << " t=" << t << " infile=" << infile << " outfile=" << outfile << endl;
        measure_graph_building_time(k, h, w, t,a, metagenomic, infile, outfolder, outfile);
    }

    return 0;
}
