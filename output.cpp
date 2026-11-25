#include "output.h"

using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::cout;

// add the deconvolution hyphen to the end of the BC or BX, wherever it is in the read
static std::string extend_barcode(std::string header, int extension)
{
    const char* whitespace = " \t\r\n\v\f";
    auto field = header.find_first_of(whitespace);

    while (field != std::string::npos) {
        field = header.find_first_not_of(whitespace, field);
        if (field == std::string::npos) {
            break;
        }

        auto end = header.find_first_of(whitespace, field);
        if (end == std::string::npos) {
            end = header.size();
        }

        if (end - field > 5 &&
            (header.compare(field, 5, "BX:Z:") == 0 ||
             header.compare(field, 5, "BC:Z:") == 0)) {
            header.insert(end, "-" + std::to_string(extension));
            break;
        }

        field = end;
    }

    return header;
}

void output(string inputFile, string outputFile, vector<Read>& reads, int min_length){

    char format = '@'; //a character to keep track of whether the input file is a fasta or a fastq
    if (inputFile.substr(inputFile.size()-6,6) == ".fasta" || inputFile.substr(inputFile.size()-3,3) == ".fa"){
        format = '>';
    }

    ifstream in(inputFile);
    if (!in){cout << "problem reading files in output_files, while trying to read " << inputFile << '\n';}

    ofstream out(outputFile);
    if (!out){
        cout << "I cannot open the output file, there may be a problem in the path of : " << outputFile << '\n';
    }


    long int sequenceID = 0;

    string tag;

    string lastnameofsequence;
    string line;
    vector<string> buffer;

    while(getline(in, line)){

        if (line[0] == format){

            //then first we append the last read we saw

            //to_deconvolve checks if the line should be deconvolved and if so returns the tag
            if (to_deconvolve(buffer, format, min_length, tag)){
                const auto key = mate_key(buffer[0], format);

                if (key == lastnameofsequence){ //if it has the same name as the previous read, it means it is paired
                    --sequenceID;
                }
                else{
                    lastnameofsequence = key;
                }

                const auto extension = reads[sequenceID].barcode_extension;
                out << extend_barcode(buffer[0], extension)
                    << '\t' << extension << '\n';
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
    if (to_deconvolve(buffer, format, min_length, tag)) {
        const auto key = mate_key(buffer[0], format);

        if (key == lastnameofsequence) {
            --sequenceID;
        } else {
            lastnameofsequence = key;
        }

        const auto extension = reads[sequenceID].barcode_extension;
        out << extend_barcode(buffer[0], extension)
            << '\t' << extension << '\n';

        ++sequenceID;
    }

}

void fasta_output(string inputFile, string outputFile, vector<Read>& reads, int min_length){

    char format = '@'; //a character to keep track of whether the input file is a fasta or a fastq
    if ((inputFile.size() > 5 && inputFile.substr(inputFile.size()-6,6) == ".fasta") || inputFile.substr(inputFile.size()-3,3) == ".fa"){
        format = '>';
    }

    ifstream in(inputFile);
    if (!in){cout << "problem reading files in output_files, while trying to read " << inputFile << '\n';}

    ofstream out(outputFile);
    if (!out){
        cout << "I cannot open the output file, there may be a problem in the path of : " << outputFile << '\n';
    }


    long int sequenceID = 0;

    string tag;

    string lastnameofsequence;
    string line;
    vector<string> buffer;

    while(getline(in, line)){

        if (line[0] == format){

            //then first we append the last read we saw

            //to_deconvolve checks if the line should be deconvolved and if so returns the tag
            if (to_deconvolve(buffer, format, min_length, tag)) {
                const auto key = mate_key(buffer[0], format);

                if (key == lastnameofsequence) {
                    --sequenceID;
                } else {
                    lastnameofsequence = key;
                }

                out << extend_barcode(buffer[0], reads[sequenceID].barcode_extension) << '\n';
                ++sequenceID;
            } else if (buffer.size()>=1){
                out << buffer[0] << '\n';
            }

            //output the rest of the lines too
            for (int i = 1 ; i < buffer.size() ; i++){
                out << buffer[i] << '\n';
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

        if (mate_key(buffer[0], format) == lastnameofsequence){ //if it has the same name as the previous read, it means it is paired
            sequenceID --;
            out << buffer[0] << "-" << reads[sequenceID].barcode_extension << '\n';
        }
        else {
            out << buffer[0] << "-" << reads[sequenceID].barcode_extension << '\n';

            lastnameofsequence = mate_key(buffer[0], format);
        }
        sequenceID++;

    }
    else if (buffer.size()>1){
        out << buffer[0] << '\n';
    }

    //output the rest of the lines too
    for (int i = 1 ; i < buffer.size() ; i++){
        out << buffer[i] << '\n';
    }

}
