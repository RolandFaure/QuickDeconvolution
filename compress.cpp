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

//returns in a deterministic fashion at least w of the reads, on average 1 / 2^hardness
//the algorithm looks if the first bit of a letter, the l
void Sequence::minimisers(int hardness, int k, int w, vector<Sequence> &minis){

    vector<bool> criterion (hardness, false);
    vector <bool> emptyWindows (this->size()-w+1, true);
    bool cont;

    do{

        int lastM = 0;
        //in empty windows, look for kmers matching the criterion
        for (int i = 0 ; i<=this->size()-k ; i++){

            int emptyWindowsSize = emptyWindows.size();
            if ((i>=emptyWindowsSize&&emptyWindows[emptyWindowsSize-1]) || (i<emptyWindowsSize&& emptyWindows[i])){

                bool good = true;
                short index = 0;
                while (good && index < hardness){
                    good = good && (s[i*2+index]==criterion[index]);
                    index ++;
                }

                if (good) { //found a minimiser !
                    minis.push_back(this->subseq(i, k));
                    for (int j = std::max(lastM, i-w) ; j <= std::min(i, emptyWindowsSize) ; j++){
                        emptyWindows[j] = false;
                    }
                    lastM = i;
                }
            }

        }

        //check if there is still an empty window
        cont = false;
        short index = 0;
        while (!cont && index<emptyWindows.size()){
            cont = cont || emptyWindows[index];
            index ++;
        }


        //update criterion by simple incrementing
        bool extra = true;
        int pos = hardness-1;
        while (extra && pos>=0){
            if (criterion[pos]){
                criterion[pos] = false;
            }
            else{
                criterion[pos] = true;
                extra = false;
            }
            pos--;
       }

    } while (cont);

}

