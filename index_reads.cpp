#include "index_reads.h"

#include <chrono>

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::pair;
using std::string;
using std::ofstream;
using std::ifstream;
using std::array;
using robin_hood::unordered_map;
using namespace std::chrono;


//the index maps a kmer to all reads beginning or ending with the kmer or its reverse complement
//the order is :
//begins by kmer : cell 0
//ends by kmer : cell 1
//ends by reverse complement : cell 3
//begins by reverse complement : cell 4

//to compress sequencing data, reads and kmers are represented by vectors of bool

void index_reads(int k, int w, string fileReads, unordered_map<Sequence, array<vector<Hit>, 4>, Sequence::HashFunction> &index, vector<pair<string, vector<long long int>>> &readClouds, vector <Read> &allreads){
	
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
    bool notag = false; //bool alerting when there is no tag attached to a read
	while(getline(in, line)){
		
        if (line[0] == '>'){
			//here looking at the name of sequence and the tag
			nameofsequence = line.erase(0,1);
			tag = get_tag(nameofsequence); //this tag is a string, as contained in a fasta: we're going to convert it into a long int, this will be much more efficient

            notag = false;
            if (tag == ""){
                notag = true;
            }
            else if (tagIDs.find(tag) == tagIDs.end()){ // if true, tag is not in the keys of tagIDs
				tagIDs[tag] = readClouds.size(); //from now on this tag ID will be associated to this tag
				tagID = readClouds.size();
                vector<long long int> s = {sequenceID};
                readClouds.push_back(std::make_pair(tag, s));

                if (readClouds.size() % 100 == 0){
                    cout << "Indexed " << readClouds.size() << " barcodes" << endl;
                }
			}
			else{
				tagID = tagIDs[tag];
                readClouds[tagID].second.push_back(sequenceID);
			}
			
			trueTag = get_true_tag(nameofsequence);

			next = true;
		}
        else if (next && !notag) {
			//here looking at the sequence itself
			next = false;
			if (k+w > line.size()){
				cout << "WARNING: the reads should be longer k+w, read " << nameofsequence << " is ignored" << endl;
                readClouds[tagID].second.pop_back();
			}
			else{
				Read r;
				r.sequence = Sequence(line);
				r.barcode = tagID;
				r.trueBarcode = trueTag;
				allreads.push_back(r);
				
				Sequence kl = r.sequence.subseq(0, k+w);
                Sequence kr = r.sequence.subseq(r.sequence.size()-k-w,k+w);
				
				pair<int, Sequence> kmerL = minimisers(kl, k, w)[0];
				
				Hit kmerL_h;
				kmerL_h.sequenceID = sequenceID;
				kmerL_h.position = kmerL.first;

                index[kmerL.second][0].push_back(kmerL_h);

				
				pair<int, Sequence> kmerR = minimisers(kr, k ,w)[0];
				Hit kmerR_h;
				kmerR_h.sequenceID = sequenceID;
				kmerR_h.position = w+k-kmerR.first;
                index[kmerR.second][1].push_back(kmerR_h);
				
				Sequence rev = kr.reverse_complement();
				pair<int, Sequence> kmerRr = minimisers(rev,k,w)[0];
				
				Hit kmerRr_h;
				kmerRr_h.sequenceID = sequenceID;
				kmerRr_h.position = kmerRr.first;
                index[kmerRr.second][2].push_back(kmerRr_h);

				
				rev = kl.reverse_complement();
				pair<int, Sequence> kmerLr = minimisers(rev,k,w)[0];
				
				Hit kmerLr_h;
				kmerLr_h.sequenceID = sequenceID;
				kmerLr_h.position = w+k-kmerLr.first;
                index[kmerLr.second][3].push_back(kmerLr_h);
			
			}
			
			sequenceID ++;
		}
	}

	auto t1 = high_resolution_clock::now();
    cout << "In index_reads, finding minimizers took me " << total_read_time/1000000 << "s out of " << duration_cast<microseconds>(t1 - t0).count()/1000000 << "s in total, to find " << index.size() << " minimizers" <<  endl;
}

//function returning all minimizers of a sequence and their positions knowing k and the windowsize


