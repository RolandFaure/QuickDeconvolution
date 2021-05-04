#include "tools.h"

#include <chrono>
using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::pair;
using std::string;
using std::ofstream;
using std::ifstream;
using robin_hood::unordered_map;
using namespace std::chrono;
//using namespace lemon;

std::string get_tag(std::string &s){
	
	string tag;
    int t =0;
    for (int i = 0; i<s.size();i++){
        /*if (s.substr(i-4,5) == "BX:Z:"){
            t = 1;
		}
        else*/ if (s[i] == ' ' or s[i] == '\t'){
            t += 1;
		}
        else if (t == 1){
			tag += s[i];
		}
	}
	return tag;
}

std::string get_true_tag(std::string &s){
	
	string tag;
	int t = 0;
	for (int i = 3; i<s.size();i++){
		if (s.substr(i-3,4) == "TBX:"){
			t = true;
		}
		else if (s[i] == ' ' or s[i] == '\\'){
			t = false;
		}
		else if (t){
			tag += s[i];
		}
	}
	return tag;
}

void export_as_SIF(std::vector<std::vector<int>> adj, std::string file){
	
	ofstream out(file);
	
	for (int i = 0 ; i<adj.size()-1 ; i++){
		for (int j = i+1 ; j<adj[0].size() ; j++){
			
			for(int k = 0 ; k<adj[i][j] ; k++){
				out<< i << "\t" << k << "\t" << j << endl;
			}
			
		}
	}
}

void export_as_CSV(std::vector<std::vector<int>> adj, std::string file){
	
	ofstream out(file);
	
	for (int i = 0 ; i<adj.size()-1 ; i++){
		for (int j = i+1 ; j<adj[0].size() ; j++){
			
			for(int k = 0 ; k<adj[i][j] ; k++){
				out<< i << "," << j << endl;
			}
			
		}
	}
}

void export_as_CSV(robin_hood::unordered_map<long int, list<int>> matching_tags, std::string file){

    ofstream out(file);

    for (robin_hood::pair<long int, list<int>> p : matching_tags){

        for (int r : p.second){
            out << std::to_string(p.first)<<"_tag,"<<r << endl;
        }

    }
}

std::string reverse_complement(std::string &s){
	
	string res = "";
	
	for (int i = s.size()-1 ; i>= 0 ; i--){
		if (s[i] == 'A'){
			res.push_back('T');
		}
		else if (s[i] == 'T'){
			res.push_back('A');
		}
		else if (s[i] == 'C'){
			res.push_back('G');
		}
		else if (s[i] == 'G'){
			res.push_back('C');
		}
	}
	return res;
}

std::vector<std::pair<int, Sequence>> minimisers(Sequence& seq, short k, short w){
	
	float total_read_time = 0;
	float total_append_time = 0;
	
	size_t seqSize = seq.size();
	//if (seqSize < k+w){cout<< "In minimisers : the input sequence is too short" << endl;}
	
	vector<int> bests(seqSize); //indices of all best kmers
	bests[0] = 0;
	int lastElement = 0; //index indicating the position of the last pertinent element 
	int firstElement = 0; //index indicating the position of the first pertinent element (which is also the best kmer at the moment) 
	
	//here indices indicate position in the sequence of nucleotide, which are position in vector<bool> seq divided by two
	vector<pair<int, Sequence>> res;
		
	//loop on all kmers
	for (int index = 1, l = seqSize+1-k; index != l; index ++){
		
		//insert the kmer among the best kmers yet
		int i = lastElement;
		bool cont = true;
		
		while (i>=firstElement && cont){
			
			bool better = seq.compare_kmers(index, bests[i], k);
			
			//cout << "Comparing " << fullnum2str(subseq(seq, index, k)) << " and " << fullnum2str(subseq(seq, bests[i],k)) << " with result " << not lexicographical_compare(seq.begin()+index*2, seq.begin()+(index+k)*2,seq.begin()+bests[i]*2, seq.begin()+(bests[i]+k)*2) << endl;
			if (better){
				
				if (bests.size() > i+1){
					bests[i+1] = index;
					lastElement = i+1;
				}
				else{
					bests.push_back(index);
					lastElement = i+1;
				}
				cont = false;
				//cout << "cont becomes false for " << fullnum2str(subseq(seq, index, k)) << endl;
			}
			i--;
		}
				
        if (index == w){ //if you have gone through the first w kmers you can choose your first minimizer
            pair<int, Sequence> p (bests[0], seq.subseq(bests[0], k));
            res.push_back(p);
        }

        // if cont == true, it means that the kmer is the better to date
        if (cont){
			
			bests[firstElement] = index;
			lastElement = firstElement;
			
            if (index>w){ //if there is already w kmers right of that, you can be sure that the new kmer is a best of a window left of that position
                pair<int, Sequence> p (index, seq.subseq(index, k));
                res.push_back(p);
            }
		}
		
        else if (index - bests[firstElement] >= w && index > w){
			while (index - bests[firstElement] >= w){
				firstElement += 1;
			}
			pair<int, Sequence> p (bests[firstElement], seq.subseq(bests[firstElement], k));
			
			res.push_back(p);
		}
		
	}
		
	//cout << "In minimizer, took me in the first phase " << total_read_time << " (among which for comparison " << total_compare_time << " and " << total_compare_time2 << ") and in the second " << total_append_time << "us out of " << duration_cast<nanoseconds>(t1 - t0).count() << "us in minimizers" << endl;

	return res;
}



//void lemon_to_csv(SmartGraph &g, std::string fileOut){
//	
//	ofstream out(fileOut);
//	int out0 = 0;
//	for(SmartGraph::EdgeIt a(g); a!=INVALID; ++a){
//		
//		out << g.id(g.u(a)) << ";" << g.id(g.v(a)) << "\n";
//	}
//}



