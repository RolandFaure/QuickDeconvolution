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
using std::set;
using robin_hood::unordered_map;
using namespace std::chrono;



//to compress sequencing data, reads and kmers are represented by vectors of bool (encapsulated in Sequence)

void index_reads(int k, int h, int w, string fileReads, vector<set<long int>> &kmers, vector<vector<long long int>> &readClouds, vector <Read> &allreads, unordered_map <string, long int> &tagIDs){
	
    long long int total_mini_time = 0;
    long long int total_read_time = 0;
	auto t0 = high_resolution_clock::now();
	
	ifstream in(fileReads);
	if (!in){cout << "problem reading files in index_reads, while trying to read " << fileReads << endl;}
	
	string nameofsequence;
	long long int sequenceID = 0; //instead of storing the name of sequence in a string we're going to store them in long long int, it will be more efficient 
	
	string tag;

	long int tagID;
	string trueTag;

    //maps all minimizing kmers to their index in kmers
    unordered_map<Sequence,long long int, Sequence::HashFunction> index;
	
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
				readClouds.push_back({sequenceID});

                if (readClouds.size() % 100 == 0){
                    cout << "Indexed " << readClouds.size() << " barcodes" << endl;
                }
			}
			else{
				tagID = tagIDs[tag];
				readClouds[tagID].push_back(sequenceID);
			}
			
			trueTag = get_true_tag(nameofsequence);

			next = true;

		}
        else if (next && !notag) {
			//here looking at the sequence itself
			next = false;
			if (k+w > line.size()){
				cout << "WARNING: the reads should be longer k+w, read " << nameofsequence << " is ignored" << endl;
				readClouds[tagID].pop_back();
			}
			else{
                auto ttt0 = high_resolution_clock::now();
				Read r;
                r.barcode = tagID;
				
                Sequence s (line);
                Sequence rev = s.reverse_complement();

                vector<Sequence> minis;
                s.minimisers(h, k, w, minis);
				
				auto ttt1 = high_resolution_clock::now();
                total_mini_time += duration_cast<nanoseconds>(ttt1 - ttt0).count();
				
                //for all minimisers in the sequence, add the tag to the kmer
                for (Sequence mini : minis){

                    auto tt0 = high_resolution_clock::now();
                    if (index.find(mini) != index.end()){
                        long int place = index[mini];
                        r.minis.push_back(place);
                        kmers[place].insert(r.barcode); //each kmer contains the set of all barcodes containing this kmer
                    }
                    else {
                        index[mini] = kmers.size();
                        r.minis.push_back(kmers.size());
                        kmers.push_back({r.barcode});
                    }
                    auto tt1 = high_resolution_clock::now();
                    total_read_time += duration_cast<nanoseconds>(tt1 - tt0).count();

                }

                minis = {};
                rev.minimisers(h, k, w, minis);
                //for all minimisers in the sequence, add the tag to the kmer
                for (Sequence mini : minis){

                    auto tt0 = high_resolution_clock::now();
                    if (index.find(mini) != index.end()){
                        long int place = index[mini];
                        r.minis.push_back(place);
                        kmers[place].insert(r.barcode); //each kmer contains the list of all barcodes containing this kmer
                    }
                    else {
                        index[mini] = kmers.size();
                        r.minis.push_back(kmers.size());
                        kmers.push_back({r.barcode});
                    }
                    auto tt1 = high_resolution_clock::now();
                    total_read_time += duration_cast<nanoseconds>(tt1 - tt0).count();

                }

                allreads.push_back(r);
			
			}
			
			sequenceID ++;
		}
	}
    auto t1 = high_resolution_clock::now();
    cout << "In index_reads, finding minimizers took me " << total_mini_time/1000000000 << "s and indexing reads " << total_read_time/1000000000 << "s out of " << duration_cast<microseconds>(t1 - t0).count()/1000000 << "s in total" <<  endl;
}

//function returning all minimizers of a sequence and their positions knowing k and the windowsize


