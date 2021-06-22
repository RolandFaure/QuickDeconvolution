#include "output.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::cout;

void output(string inputFile, string outputFile, vector<Read>& reads){

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

    string line;
    bool notag = false; //bool alerting when there is no tag attached to a read
    while(getline(in, line)){

        if (line[0] == format){

            //here looking at the name of sequence and the tag
            nameofsequence = line.erase(0,1);

            out << nameofsequence << "\t" << reads[sequenceID].barcode_extension << endl;

            sequenceID++;
        }

    }

}
