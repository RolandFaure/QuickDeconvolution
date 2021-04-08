#include "index_reads.h"

#include <chrono>

using namespace std;
using namespace std::chrono;

//the index maps a kmer to all reads beginning or ending with the kmer or its reverse complement
//the order is :
//begins by kmer : cell 0
//ends by kmer : cell 1
//ends by reverse complement : cell 3
//begins by reverse complement : cell 4

//to compress sequencing data, reads and kmers are represented by vectors of bool

void index_reads(int k, int w, string fileReads, unordered_map<Sequence, vector<vector<Hit>>, Sequence::HashFunction> &index, vector<vector<long long int>> &readClouds, vector <Read> &allreads){
	
	float total_read_time = 0;
	auto t0 = high_resolution_clock::now();
	
	ifstream in(fileReads);
	if (!in){cout << "problem reading files in index_reads, while trying to read " << fileReads << endl;}
	
	string nameofsequence;
	long long int sequenceID = 0; //instead of storing the name of sequence in a string we're going to store them in long long int, it will be more efficient 
	
	string tag;
	unordered_map <string, long int> tagIDs;
	long int tagID;
	string trueTag;
	
	string line;
	bool next = false;
	while(getline(in, line)){
		
		if (line[0] == '@'){
			//here looking at the name of sequence and the tag
			nameofsequence = line.erase(0,1);
			tag = get_tag(nameofsequence); //this tag is a string, as contained in a fasta: we're going to convert it into a long int, this will be much more efficient
			
			if (tagIDs.find(tag) == tagIDs.end()){ // if true, tag is not in the keys of tagIDs
				tagIDs[tag] = readClouds.size(); //from now on this tag ID will be associated to this tag
				tagID = readClouds.size();
				readClouds.push_back({sequenceID});
			}
			else{
				tagID = tagIDs[tag];
				readClouds[tagID].push_back(sequenceID);
			}
			
			trueTag = get_true_tag(nameofsequence);

			next = true;
		}
		else if (next) {
			//here looking at the sequence itself
			next = false;
			if (k+w > line.size()){
				cout << "WARNING: the reads should be longer k+w, read " << nameofsequence << " is ignored" << endl;
				readClouds[tagID].pop_back();
			}
			else{
				Read r;
				r.sequence = Sequence(line);
				r.barcode = tagID;
				r.trueBarcode = trueTag;
				allreads.push_back(r);
				
				Sequence kl = r.sequence.subseq(0, k+w);
				Sequence kr = r.sequence.subseq(r.sequence.size()/2-k-w,k+w);
				
				auto ttt0 = high_resolution_clock::now();
				pair<int, Sequence> kmerL = minimisers(kl, k, w)[0];
				
				auto ttt1 = high_resolution_clock::now();
				total_read_time += duration_cast<microseconds>(ttt1 - ttt0).count();
				
				Hit kmerL_h;
				kmerL_h.sequenceID = sequenceID;
				kmerL_h.position = kmerL.first;
				if (index.count(kmerL.second) > 0) {
					index[kmerL.second][0].push_back(kmerL_h);
				}
				else{
					index[kmerL.second] = {{kmerL_h},{},{},{}};
				}
				
				ttt0 = high_resolution_clock::now();
				pair<int, Sequence> kmerR = minimisers(kr, k ,w)[0];
				ttt1 = high_resolution_clock::now();
				total_read_time += duration_cast<microseconds>(ttt1 - ttt0).count();
				
				Hit kmerR_h;
				kmerR_h.sequenceID = sequenceID;
				kmerR_h.position = w+k-kmerR.first;
				if (index.count(kmerR.second) > 0) {
					index[kmerR.second][1].push_back(kmerR_h);
				}
				else{
					index[kmerR.second] = {{},{kmerR_h},{},{}};
				}
				
				ttt0 = high_resolution_clock::now();
				Sequence rev = kr.reverse_complement();
				pair<int, Sequence> kmerRr = minimisers(rev,k,w)[0];
				ttt1 = high_resolution_clock::now();
				total_read_time += duration_cast<microseconds>(ttt1 - ttt0).count();
				
				Hit kmerRr_h;
				kmerRr_h.sequenceID = sequenceID;
				kmerRr_h.position = kmerRr.first;
				if (index.count(kmerRr.second) > 0) {
					index[kmerRr.second][2].push_back(kmerRr_h);
				}
				else{
					index[kmerRr.second] = {{},{},{kmerRr_h},{}};
				}
				
				ttt0 = high_resolution_clock::now();
				rev = kl.reverse_complement();
				pair<int, Sequence> kmerLr = minimisers(rev,k,w)[0];
				ttt1 = high_resolution_clock::now();
				total_read_time += duration_cast<microseconds>(ttt1 - ttt0).count();
				
				Hit kmerLr_h;
				kmerLr_h.sequenceID = sequenceID;
				kmerLr_h.position = w+k-kmerLr.first;
				if (index.count(kmerLr.second) > 0) {
					index[kmerLr.second][3].push_back(kmerLr_h);
				}
				else{
					index[kmerLr.second] = {{},{},{},{kmerLr_h}};
				}
			
			}
			
			sequenceID ++;
		}
	}
	auto t1 = high_resolution_clock::now();
	//cout << "In index_reads, finding minimizers took me " << total_read_time/1000000 << "s out of " << duration_cast<microseconds>(t1 - t0).count()/1000000 << "s in total" <<  endl;
}

//function returning all minimizers of a sequence and their positions knowing k and the windowsize

