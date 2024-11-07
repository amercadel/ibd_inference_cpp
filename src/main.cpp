#include <iostream>
#include "vcf.hpp"
#include "RPBWT.hpp"
#include "PBWT.hpp"

int main(int argc, char** argv){
    std::string input_vcf = std::string(argv[1]);
    std::string output_vcf = std::string(argv[2]);
    int checkpoint = std::atoi(argv[3]);
    std::string genetic_map = std::string(argv[4]);
    int L = 20;
    int W = 50;
    int G = 1;
    double rho = 0.05;
    int ret = runRPBWT(input_vcf, output_vcf, checkpoint);
    ret = runPBWT(input_vcf, output_vcf, genetic_map, checkpoint, L, W, G, rho);
}