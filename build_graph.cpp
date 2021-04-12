#include "build_graph.h"
#include "cluster_graph.h"

//#include <lemon/smart_graph.h> //graph library
//#include <lemon/list_graph.h>

using namespace std;
using namespace std::chrono; 
//using namespace lemon;

//the function takes as an input the list of all reads having the same tag
vector<int> build_graph(int k, int w, long int tagCloud, std::vector<long long int>& readCloud, std::vector <Read> &reads, std::unordered_map<Sequence, vector<std::vector<Hit>>, Sequence::HashFunction> &index){
	
	long int mini_time = 0;
	auto t0 = high_resolution_clock::now();
		
	unordered_map<long int, list<int>> matching_tags; //maps to each tagID the index of the reads of the cloud aligning with this tag
	
	float total_read_time = 0;
	
	//int r = 0;
	for(int r = 0, sizer = readCloud.size(); r<sizer ; r++){

		long long int name = readCloud[r];
		
		Sequence read = reads[name].sequence; //name is the name of the read, read is the sequence
		//cout << "\nLooking at read " << fullnum2str(read) << endl;
		
		auto tt0 = high_resolution_clock::now();
		std::vector<std::pair<int, Sequence>> mini = minimisers(read, k, w); //mini is the list of all minimisers of the read
		auto tt1 = high_resolution_clock::now();
		mini_time += duration_cast<nanoseconds>(tt1-tt0).count();
		
		for (auto m : mini){
			
			int s = m.first; //postion of the minimizer
			Sequence sub = m.second; //sequence of the minimizing kmer
			
			auto ttt0 = high_resolution_clock::now();
			auto candidates = index[sub];
			auto ttt1 = high_resolution_clock::now();
			total_read_time += duration_cast<microseconds>(ttt1 - ttt0).count();

			
			//cout << "Minimizer " << fullnum2str(m.second) << endl;
			
			//if candidates.size()!=4, it means there are no corresponding kmers
			if (candidates.size() == 4){
				
				Hit candidateL;
				for (int c=0, sizec = candidates[0].size(); c<sizec; c++){
					
					candidateL = candidates[0][c];
					
					Read candid = reads[candidateL.sequenceID];
					Sequence candidate = candid.sequence;
					long int tag = candid.barcode;
					int alignmentLength = read.size()-s+candidateL.position;
					

					if (candid.barcode != tagCloud && alignmentLength<=read.size() && read.subseq(s-candidateL.position, alignmentLength) == candidate.subseq(0 ,alignmentLength)){
						
						//cout << "Match left with " << fullnum2str(candid.sequence) << endl;						
						matching_tags[tag].push_back(r);
					}
				}
				
				
				Hit candidateR;
				for (int c=0, sizec = candidates[1].size(); c<sizec; c++){
					
					candidateR = candidates[1][c];
					
					Read candid = reads[candidateR.sequenceID];
					Sequence candidate = candid.sequence;
					long int tag = candid.barcode;					
					int alignmentLength = min(s+candidateR.position, int(candidate.size()));

					
					//cout << "Candidate right with " << fullnum2str(candidate) << ", comparing " << fullnum2str(subseq(candidate, candidate.size()/2-alignmentLength, alignmentLength)) << " and " << fullnum2str(subseq(read, s+candidateR.position-alignmentLength,alignmentLength)) << " (alignment length: " << alignmentLength << "), minimizer: " << fullnum2str(sub)  << endl;

					if (candid.barcode != tagCloud && s+candidateR.position<=candidate.size() && candidate.subseq(candidate.size()-alignmentLength, alignmentLength) == read.subseq(s+candidateR.position-alignmentLength,alignmentLength)){
						
						//cout << "Match right with " << fullnum2str(candid.sequence) << endl;						
						matching_tags[tag].push_back(r);
					}
				}

				Hit candidateRr;
				for (int c=0, sizec = candidates[2].size(); c<sizec; c++){
					
					candidateRr = candidates[2][c];
					Read candid = reads[candidateRr.sequenceID];
					Sequence candidate = candid.sequence.reverse_complement();
					long int tag = candid.barcode;
					
					int alignmentLength = read.size()-s+candidateRr.position;
					//cout << "candidateRr : " << fullnum2str(candidate) << " (" << fullnum2str(candid.sequence) << "), comparing " << fullnum2str(subseq(candidate, 0, alignmentLength)) << " and " << fullnum2str(subseq(read,s-candidateRr.position,alignmentLength)) << endl;

					if (candid.barcode != tagCloud && alignmentLength<=candidate.size() && candidate.subseq(0, alignmentLength) == read.subseq(s-candidateRr.position,alignmentLength)){
						
						//cout << "Match !" << endl;						
						matching_tags[tag].push_back(r);
					}
				}

				Hit candidateLr;
				for (int c=0, sizec = candidates[3].size(); c<sizec; c++){
					
					candidateLr = candidates[3][c];
					Read candid = reads[candidateLr.sequenceID];
					Sequence candidate = candid.sequence.reverse_complement(); // we have a sequence which ends with the kmer
					long int tag = candid.barcode;
					
					int alignmentLength = min(s+candidateLr.position, int(candidate.size()));
					
					if (candid.barcode != tagCloud && s+candidateLr.position<=read.size() && candidate.subseq(candidate.size()-alignmentLength, alignmentLength) == read.subseq(s+candidateLr.position-alignmentLength,alignmentLength)){
						
						//cout << "Match !" << endl;													
						matching_tags[tag].push_back(r);
					}
				}
			}
            else { // in this case, by trying to access index[sub] above we've created a useless key sub in index, which will take much space in memory
                index.erase(sub);
            }
		}
		//r++;
	}

	auto t1 = high_resolution_clock::now();
	
	vector<int> clusters (readCloud.size(), -1);
	cluster_graph(matching_tags, clusters);
	
//	int n = 0;
////	cout << "sequence of read 0 : " << fullnum2str(reads[0].sequence) << endl;
////	cout << "sequence of read 0' : " << fullnum2str(reads[1].sequence) << endl;
//	if (tagCloud == 0){
//		for (pair<long int, list<int>> matchs : matching_tags){
//			for (int i: matchs.second){
//				if (readCloud[i] == 0){
//					cout << "Read 0 matched with tag " << matchs.first << endl;
//					n++;
//				}
//				if (readCloud[i] == 1){
//					cout << "Read 0' matched with tag " << matchs.first << endl;
//					n++;
//				}
//			}
//		}
//	}
//	cout << "reads 0 matched with " << n << " tags" << endl;

	auto t2 = high_resolution_clock::now();
	
	//cout << "push_back time : " << total_read_time << endl;
	//cout << "Alignement time : " << duration_cast<nanoseconds>(t1 - t0).count() <<"us, total read accession : " << total_read_time << ", total minimizer time : " << mini_time << endl;
	return clusters;
}
