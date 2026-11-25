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
#include <thread>
#include <mutex>
#include <atomic> //for multi-threading
#include <fstream>
#include <stdexcept>
#include <string>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "robin_hood.h"

#include "compress.h"

struct Hit{
	long long int sequenceID;
	int position;
};

//struct Read{
//    std::vector<std::vector<long int>> minis; //one vector for each thread (in the simplest case, juste a vector)
//    long int barcode;
//};

inline bool is_gzip_path(const std::string& path)
{
    const auto n = path.size();
    return n >= 3 && path[n - 3] == '.' &&
           (path[n - 2] == 'g' || path[n - 2] == 'G') &&
           (path[n - 1] == 'z' || path[n - 1] == 'Z');
}
inline char sequence_format(const std::string& path)
{
    const auto end = path.size() - (is_gzip_path(path) ? 3 : 0);
    const bool fasta =
        (end >= 6 && path.compare(end - 6, 6, ".fasta") == 0) ||
        (end >= 3 && path.compare(end - 3, 3, ".fa") == 0);
    return fasta ? '>' : '@';
}

bool to_deconvolve(std::vector<std::string> &buffer, char format, int min_length, std::string &tag);

std::string get_tag(std::string &s, char format);
std::string get_true_tag(std::string &s);
void export_as_SIF(std::vector<std::vector<int>> adj, std::string file);
void export_as_CSV(std::vector<std::vector<int>> &adj, std::string file,  std::string fileNode, std::vector<int> &clusters);
void export_as_CSV(robin_hood::unordered_map<long int, std::unordered_set<int>> matching_tags, std::string file);
std::string reverse_complement(std::string &s);
std::string mate_key(const std::string& header, char format);

std::ostream& operator<< (std::ostream& out, const std::vector<int>& v);
void minimisers(Sequence& seq, short k, short w, std::vector<Sequence> &res);


inline void open_sequence_input(
    const std::string& path,
    std::ifstream& file,
    boost::iostreams::filtering_istream& in)
{
    file.open(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open input file " + path);
    }
    file.exceptions(std::ios::badbit);

    if (is_gzip_path(path)) {
        in.push(boost::iostreams::gzip_decompressor());
    }
    in.push(file);
    in.exceptions(std::ios::badbit);
}

inline void open_sequence_output(
    const std::string& path,
    std::ofstream& file,
    boost::iostreams::filtering_ostream& out)
{
    file.open(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open output file " + path);
    }

    if (is_gzip_path(path)) {
        out.push(boost::iostreams::gzip_compressor());
    }
    out.push(file);
    out.exceptions(std::ios::badbit | std::ios::failbit);
}

inline void finish_sequence_output(
    const std::string& path,
    std::ofstream& file,
    boost::iostreams::filtering_ostream& out)
{
    out.reset();  // Finish the filter chain, including the gzip trailer.
    file.close();
    if (!file) {
        throw std::runtime_error("Could not finish writing " + path);
    }
}

//std::vector<bool> subseq(std::vector<bool> &seq, int start, int length);

//void lemon_to_csv(lemon::SmartGraph &g, std::string fileOut);

#endif
