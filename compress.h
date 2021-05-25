#ifndef COMPRESS_H
#define COMPRESS_H

#include <string>
#include <vector>
#include <iostream>
 
class Sequence{
	
public :
	
	Sequence();
	Sequence(std::string &inputSequence);
	Sequence(std::vector<bool> &inputVector);
	
	Sequence reverse_complement() const; //returns a reverse complement sequence
	Sequence subseq(int start, int length);
	std::string str() const; //returns a string of ACGT
	size_t size();
		
    void minimisers(int hardness, int k, int w, std::vector<std::vector<Sequence>> &minis);
	bool compare_kmers(int start1, int start2, int k); //specific function to compare leicographically to kmer within a sequence (useful to find minimizers)
	
	struct HashFunction
	{
		size_t operator()(const Sequence& seq) const
		{
			std::hash<std::vector<bool>> h;
			return h(seq.s);
		}
	};
    size_t hash();
		
private :
	
	std::vector<bool> s;
	
//    friend std::ostream& operator<< (std::ostream& stream, Sequence const& sequence);
	friend bool operator==(Sequence const &a , Sequence const &b);

}; 


//std::vector<bool> fullstr2num(const std::string& str);
std::string fullnum2str(std::vector<bool> num);
//
//std::vector<bool> reverse_complement(std::vector<bool> &seq);

#endif
