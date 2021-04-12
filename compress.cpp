#include "compress.h"

using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::string;

Sequence::Sequence(){
	
}

Sequence::Sequence(string& inputSequence){

	for(uint i(0);i<inputSequence.size();i++){
		switch (inputSequence[i]){
		  case 'A': s.push_back(false); s.push_back(false);break;
		  case 'C': s.push_back(false); s.push_back(true);break;
		  case 'G': s.push_back(true); s.push_back(false);break;
		  default: s.push_back(true); s.push_back(true);break;
		}
	}
}

Sequence::Sequence(vector<bool> &inputVector){
	
	s = inputVector;
}

string Sequence::str() const{
	
	string str(s.size()/2, 'N');
	uint j = 0;
	for(uint i(0);i<s.size();i+=2){
		if(s[i]){
		  if(s[i+1]){
			  str[j] = 'T';
		  }else{
			str[j] = 'G';
		  }
		}else{
		  if(s[i+1]){
			str[j] = 'C';
		  }else{
			str[j] = 'A';
		  }
		}
		j++;
  }
  return str;
	
}

Sequence Sequence::reverse_complement() const{
	
	std::vector<bool> res;
	
	for (int i = s.size()-1 ; i>0 ; i-=2){
		res.push_back(not s[i-1]);
		res.push_back(not s[i]);
	}
	
	Sequence r (res);
	return r;
}

//takes a subset of a sequence, with argument the position on the sequence and the number of nucleotides
Sequence Sequence::subseq(int start, int length){
	
//this verification might be useful for debugging, but it slows down the program A LOT
//	if (start*2+length*2 > seq.size()){
//		cout << "ERROR in subseq, the sequence is too short " << seq.size() << " " << start << " " << length << endl;
//	}
	vector<bool>::const_iterator st = s.begin() + 2*start;
	vector<bool>::const_iterator f = s.begin() + 2*length+2*start;
	vector<bool> res(st, f);
	
	//auto s = seq.begin();
	return Sequence(res);
}

bool operator==(Sequence const &a , Sequence const &b){
	return a.s == b.s;
}


size_t Sequence::size(){
	return size_t(s.size()/2);
}

//returns true if the range [start1,start1+k) less than [start2, start2+k)
bool Sequence::compare_kmers(int start1, int start2, int k){
	return not lexicographical_compare(s.begin()+start1*2, s.begin()+(start1+k)*2,s.begin()+start2*2, s.begin()+(start2+k)*2);
}
//std::vector<bool> fullstr2num(const string& str) {
//  std::vector<bool> res;
//  for(uint i(0);i<str.size();i++){
//    switch (str[i]){
//      case 'A':res.push_back(false);res.push_back(false);break;
//      case 'C':res.push_back(false);res.push_back(true);break;
//      case 'G':res.push_back(true);res.push_back(false);break;
//      default:res.push_back(true);res.push_back(true);break;
//    }
//  }
//  return res;
//}

//std::string fullnum2str(vector<bool> num) {
//  string str(num.size()/2, 'N');
//  uint j = 0;
//  for(uint i(0);i<num.size();i+=2){
//    if(num[i]){
//      if(num[i+1]){
//          str[j] = 'T';
//      }else{
//        str[j] = 'G';
//      }
//    }else{
//      if(num[i+1]){
//        str[j] = 'C';
//      }else{
//        str[j] = 'A';
//      }
//    }
//    j++;
//  }
//  return str;
//}

//std::vector<bool> reverse_complement(std::vector<bool> &seq){
//	
//	std::vector<bool> res;
//	seq.flip();
//	
//	for (int i = seq.size()-1 ; i>0 ; i-=2){
//		res.push_back(seq[i-1]);
//		res.push_back(seq[i]);
//	}
//	return res;
//}
