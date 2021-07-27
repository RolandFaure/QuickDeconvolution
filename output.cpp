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


    long int sequenceID = 0;

    string tag;

    string lastnameofsequence;
    string line;
    vector<string> buffer;

    while(getline(in, line)){

        if (line[0] == format){

            //then first we append the last read we saw

            if (to_deconvolve(buffer, format, min_length, tag)){ //to_deconvolve checks if the line should be deconvolved and if so returns the tag

                if (buffer[0] == lastnameofsequence){ //if it has the same name as the previous read, it means it is paired
                    sequenceID --;
                    out << buffer[0] << "\t" << reads[sequenceID].barcode_extension << endl;
                }
                else{
                    out << buffer[0] << "\t" << reads[sequenceID].barcode_extension << endl;

                    lastnameofsequence = buffer[0];
                }
                sequenceID++;

            }

            //then we reset the buffer
            buffer = {line};

        }
        else {
            buffer.push_back(line);
        }

    }

    //parse the last read
    if (to_deconvolve(buffer, format, min_length, tag)){

        if (buffer[0] == lastnameofsequence){ //if it has the same name as the previous read, it means it is paired
            out << buffer[0] << "\t" << reads[sequenceID].barcode_extension << endl;
        }
        else{
            out << buffer[0] << "\t" << reads[sequenceID].barcode_extension << endl;

            sequenceID ++;

            lastnameofsequence = buffer[0];
        }


    }

}

void fasta_output(string inputFile, string outputFile, vector<Read>& reads, int min_length){

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


    long int sequenceID = 0;

    string tag;

    string lastnameofsequence;
    string line;
    vector<string> buffer;

    while(getline(in, line)){

        if (line[0] == format){

            //then first we append the last read we saw

            if (to_deconvolve(buffer, format, min_length, tag)){ //to_deconvolve checks if the line should be deconvolved and if so returns the tag

                if (buffer[0] == lastnameofsequence){ //if it has the same name as the previous read, it means it is paired
                    sequenceID --;
                    out << buffer[0] << "-" << reads[sequenceID].barcode_extension << endl;
                }
                else{
                    out << buffer[0] << "-" << reads[sequenceID].barcode_extension << endl;

                    lastnameofsequence = buffer[0];
                }
                sequenceID++;

            }
            else if (buffer.size()>1){
                out << buffer[0] << endl;
            }

            //output the rest of the lines too
            for (int i = 1 ; i < buffer.size() ; i++){
                out << buffer[i] << endl;
            }

            //then we reset the buffer
            buffer = {line};

        }
        else {
            buffer.push_back(line);
        }

    }

    //parse the last read
    if (to_deconvolve(buffer, format, min_length, tag)){ //to_deconvolve checks if the line should be deconvolved and if so returns the tag

        if (buffer[0] == lastnameofsequence){ //if it has the same name as the previous read, it means it is paired
            sequenceID --;
            out << buffer[0] << "-" << reads[sequenceID].barcode_extension << endl;
        }
        else {
            out << buffer[0] << "-" << reads[sequenceID].barcode_extension << endl;

            lastnameofsequence = buffer[0];
        }
        sequenceID++;

    }
    else if (buffer.size()>1){
        out << buffer[0] << endl;
    }

    //output the rest of the lines too
    for (int i = 1 ; i < buffer.size() ; i++){
        out << buffer[i] << endl;
    }

}
