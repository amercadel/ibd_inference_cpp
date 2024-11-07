#pragma once
#include<fstream>
#include <string>

void skipMeta(std::ifstream& in);
void getM(std::ifstream& in);
int runRPBWT(std::string input, std::string output, int checkpoint);