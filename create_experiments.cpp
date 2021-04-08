#include "create_experiments.h"
#include "generate_random_sequence.h"
#include "draw_10x.h"

using namespace std;
using namespace chrono;

void create_exps(){
	
	char alphabet[4] = {'A', 'C', 'G', 'T'};
	
//	chromosomes_to_file(100000, alphabet, 4, 1, "eval/genome_100kb.fasta");
	chromosomes_to_file(1000000, alphabet, 4, 1, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_1Mb.fasta");
	chromosomes_to_file(10000000, alphabet, 4, 1, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_10Mb.fasta");
	chromosomes_to_file(100000000, alphabet, 4, 1, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_100Mb.fasta");
//	chromosomes_to_file(100000000, alphabet, 4, 10, "eval/genome_1Gb.fasta");
//	
//	draw_fragments(70, 150, 25, 0.2, 9000, 11000, 4, "eval/genome_100kb.fasta", "eval/reads_100kb_cov25_redundance4.fasta");
//	draw_fragments(15, 40, 50, 0.2, 2500, 7500, 4, "eval/genome_100kb.fasta", "eval/reads_100kb_cov50_redundance4.fasta");
//	draw_fragments(15, 40, 50, 0.2, 2500, 7500, 15, "eval/genome_100kb.fasta", "eval/reads_100kb_cov25_redundance15.fasta");
//	
	draw_fragments(150, 400, 25, 0.2, 25000, 75000, 4, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_1Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_1Mb_cov25_redundance4.fasta");
	draw_fragments(150, 400, 50, 0.2, 25000, 75000, 4, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_1Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_1Mb_cov50_redundance4.fasta");
	draw_fragments(150, 400, 50, 0.2, 25000, 75000, 15, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_1Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_1Mb_cov50_redundance15.fasta");
	
	draw_fragments(150, 400, 25, 0.2, 25000, 75000, 4, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_10Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_10Mb_cov25_redundance4.fasta");
	draw_fragments(150, 400, 50, 0.2, 25000, 75000, 4, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_10Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_10Mb_cov50_redundance4.fasta");
	draw_fragments(150, 400, 50, 0.2, 25000, 75000, 15, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_10Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_10Mb_cov50_redundance15.fasta");
	
	draw_fragments(150, 400, 25, 0.2, 25000, 75000, 4, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_100Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_100Mb_cov25_redundance4.fasta");
	draw_fragments(150, 400, 50, 0.2, 25000, 75000, 4, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_100Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_100Mb_cov50_redundance4.fasta");
	draw_fragments(150, 400, 50, 0.2, 25000, 75000, 15, "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/genome_100Mb.fasta", "/home/zaltabar/Documents/Ecole/X/4A/stage_M2/code/10xSplitter/eval/reads_100Mb_cov50_redundance15.fasta");
	
//	draw_fragments(150, 400, 25, 0.2, 25000, 75000, 4, "eval/genome_1Gb.fasta", "eval/reads_1Gb_cov25_redundance4.fasta");
//	draw_fragments(150, 400, 50, 0.2, 25000, 75000, 4, "eval/genome_1Gb.fasta", "eval/reads_1Gb_cov50_redundance4.fasta");
//	draw_fragments(150, 400, 50, 0.2, 25000, 75000, 15, "eval/genome_1Gb.fasta", "eval/reads_1Gb_cov50_redundance15.fasta");
}