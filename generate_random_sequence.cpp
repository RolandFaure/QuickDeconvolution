#include "generate_random_sequence.h"

#include <string>
#include <iostream>
#include <ctime> // pour tirer des nombres aléatoires
#include <cstdlib> // pour tirer des nombres aléatoires
#include <vector>
#include <fstream>

using std::cout;
using std::endl;
using std::vector;
using std::pair;
using std::string;
using std::ofstream;
using std::ifstream;

string uniform_random_sequence(int length, char alphabet[], int sizeAlphabet){

    string seq = "";
    int nb;
    srand(time(0));

    for(int i = 0 ; i < length ; i++){

        nb = rand() % sizeAlphabet;

        seq += alphabet[nb];

    }

    return seq;
}

void chromosomes_to_file(int length, char alphabet[], int sizeAlphabet, int nbChr, string fileOut){

    ofstream out(fileOut.c_str());
    if (!out){ cout << "Error while opening the file in chromosome_to_file" << endl;}

    for (int i = 0 ; i < nbChr ; i++){

        string s;
        s = uniform_random_sequence(length, alphabet, sizeAlphabet);
        out << ">chr"<<i<<"\n"<<s<<endl;
    }



}
