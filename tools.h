#ifndef TOOLS_D
#define TOOLS_D

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <list>
//#include <lemon/smart_graph.h> //graph library

#include "compress.h"

struct Hit{
	long long int sequenceID;
	int position;
};

struct Read{
	Sequence sequence;
	long int barcode;
	std::string trueBarcode;
};

std::string get_tag(std::string &s);
std::string get_true_tag(std::string &s);
void export_as_SIF(std::vector<std::vector<int>> adj, std::string file);
void export_as_CSV(std::vector<std::vector<int>> adj, std::string file);
std::string reverse_complement(std::string &s);

std::vector<std::pair<int, Sequence>> minimisers(Sequence &seq, short k, short w);

//std::vector<bool> subseq(std::vector<bool> &seq, int start, int length);

//void lemon_to_csv(lemon::SmartGraph &g, std::string fileOut);

#endif