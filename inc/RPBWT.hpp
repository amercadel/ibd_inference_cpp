#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <numeric>
#include <cassert>

using namespace std;



// skip meta-info lines and header line
void skipMeta(ifstream& in);

// find M using vcf's info field - resets input stream pointer when finished
int getM(ifstream& in);

int runRPBWT(string input_vcf, string write_to, int checkpoint);