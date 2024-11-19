#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "utils.hpp"



struct rateMapData{
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    std::vector<int> all_sites;
    std::vector<float> interpolated_cm;
    std::vector<float> interpolateVector(std::vector<int> &sites, std::vector<int> &bp_vec, std::vector<float> &cm_vec);
    void write_interpolated_map(std::string output_file);
    std::vector<float> interpolateVector(std::vector<int> &sites);
    float interpolateBasePairToGenPos(int site);
    float genPos(int site);
    int last_bp;
    float last_cm;
    int min_tree_subsample;

    void display(int idx){
        std::cout << "bp value: " << bp_vec[idx] << std::endl << "cm value: " << cm_vec[idx] << std::endl;
    };

};
rateMapData readRateMap(char *filename);
rateMapData readRateMap(char* filename, std::vector<int> &sites);

// gets genetic position based on interpolated genetic distance
// OUTPUT: a float representing a site's genetic position
float getGeneticPosition(std::vector<float> &interpolated_cm, int site_idx);
