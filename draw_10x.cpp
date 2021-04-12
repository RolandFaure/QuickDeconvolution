#include "draw_10x.h"

#include <ctime> // for generating random numbers
#include <cstdlib> // for generating random numbers
#include <chrono>

using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::string;
using std::ifstream;
using std::ofstream;

void draw_fragments(int const read_size, int insert_size, int gen_coverage, float frag_coverage, int frag_size_min, int frag_size_max, float redundance, std::string fileGenome, std::string fileReads){
	
	ifstream in(fileGenome.c_str());
	if (not in){
		cout << "Problem while reading file in draww_fragments, check file name" << endl;
	}
	
	ofstream out(fileReads);
	
	string line;
	vector<string> chromosomes;
	
	int genome_size = 0;
	
	int delta = insert_size - 2*read_size;
	
	while(getline(in, line)){
		
		if (line[0] != '>'){
			chromosomes.push_back(line);
			genome_size += line.size();
		}
	}
	

	long long int COV = 0;
	std::random_device rd;
	srand(time(0));
	
	long int nb_tags = long(genome_size * gen_coverage / (frag_coverage*0.5*(frag_size_max+frag_size_min)) / redundance);
	
	int ifragment = 0; //an index unique to each fragment
	int iread = 0; //an index unique to each read, to build the fasta
	
	//cout << genome_size << " " << gen_coverage << genome_size * gen_coverage <<endl;
	
	while (COV < genome_size / 100 * gen_coverage){
	
		// select ramdomly a chromosome
		int num_chr = rand()%(chromosomes.size());
		// extract a fragment
		int length = rand()%(int(frag_size_max-frag_size_min))+frag_size_min;
		long start = rd()%(chromosomes[num_chr].size()-length);
		
		COV += int(length*frag_coverage/100);
		
		//now extract the reads from the fragment
		//affect randomly a tag to that frag
		long num_tag = rd()%nb_tags;
		// compute the number of paired read for that fragment
		
		int nb_paired_reads = int((length*frag_coverage)/(2*read_size));
		int index_pr = 0;
		while (index_pr < nb_paired_reads){
			int gap = rand()%delta;
			int i = rand()%(length-2*read_size-gap);
			int j = i + read_size + gap;
						
			string quality = "";
			for (int i = 0;i<read_size;i++){quality += 'A';} //the line indicating the quality of the read in fastq format
			out << "@read"<< iread << "_TBX:" << ifragment << "\tBX:Z:" << num_tag << "\n" << chromosomes[num_chr].substr(start+i, read_size) << "\n+\n" << quality << endl;
			out << "@read"<< iread << "_TBX:" << ifragment<<  "\tBX:Z:" << num_tag << "\n" << chromosomes[num_chr].substr(start+j, read_size) << "\n+\n" << quality << endl;
			
			index_pr += 1;
			iread += 1;
		}
		ifragment += 1;
		
		cout << "Generated " << float(COV)/(genome_size/100 * gen_coverage) *100 << "% of fake paired 10x reads" << endl;
	}
}
