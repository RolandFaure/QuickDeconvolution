#ifndef TOOLS_D
#define TOOLS_D

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <list>
#include <set>
#include <unordered_set>
#include <chrono>

#include "robin_hood.h"

#include "compress.h"

struct Hit{
	long long int sequenceID;
	int position;
};

struct Read{
    std::vector<long int> minis;
	long int barcode;
};

//class Filtering_set{ //an unordered set, except that you add an element only if you already have seen it n-1 times

//public :

//    Filtering_set(n);
//    std::unordered_set set;

//private :

//    int n;

//};

std::string get_tag(std::string &s, char format);
std::string get_true_tag(std::string &s);
void export_as_SIF(std::vector<std::vector<int>> adj, std::string file);
void export_as_CSV(std::vector<std::vector<int>> &adj, std::string file,  std::string fileNode, std::vector<int> &clusters);
void export_as_CSV(robin_hood::unordered_map<long int, std::unordered_set<int>> matching_tags, std::string file);
std::string reverse_complement(std::string &s);

std::ostream& operator<< (std::ostream& out, const std::vector<int>& v);
void minimisers(Sequence& seq, short k, short w, std::vector<Sequence> &res);

//std::vector<bool> subseq(std::vector<bool> &seq, int start, int length);

//void lemon_to_csv(lemon::SmartGraph &g, std::string fileOut);

#endif
