#include "output.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::cout;

void output(string inputFile, string outputFile, vector<Read>& reads, int min_length){

    char format = '@'; //a character to keep track of whether the input file is a fasta or a fastq
    if (inputFile.substr(inputFile.size()-6,6) == ".fasta" || inputFile.substr(inputFile.size()-3,3) == ".fa"){
        format = '>';
    }

    ifstream in(inputFile);
    if (!in){cout << "problem reading files in output_files, while trying to read " << inputFile << endl;}

    ofstream out(outputFile);
    if (!out){
        cout << "I cannot open the output file, there may be a problem in the path of : " << outputFile << endl;
    }

    string nameofsequence;

    long int sequenceID;

    vector<string> buffer;
    string tag;

    string line;
    while(getline(in, line)){


        if (line[0] == format){

            //here looking at the name of sequence and the tag

            if (to_deconvolve(buffer, format, min_length, tag)){
                nameofsequence = line.erase(0,1);

                out << nameofsequence << "\t" << reads[sequenceID].barcode_extension << endl;

                sequenceID++;
            }

        }
        else{
            buffer.push_back(line);
        }

    }

}
