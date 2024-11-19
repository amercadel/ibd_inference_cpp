#pragma once
#include <iostream>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <string>
#include "read_rate_map.hpp"
#include "utils.hpp"
 



// The Match class represents an IBS seed. By definition it must have greater than min_seed length and more than min_markers sites within it
// It has a constructor from a string and by directly inputting its data; only the from string version is used
class Match{
    public:
        int hap1;
        int hap2;
        int start_site;
        int end_site;
        int n_sites;
        float len_cm = 0.0;
        
        Match(std::string &input_str){
            std::vector<std::string> split_str = split(input_str, '\t');
            this->hap1 = std::min(stoi(split_str[1]), stoi(split_str[2]));
            this->hap2 = std::max(stoi(split_str[1]), stoi(split_str[2]));
            this->start_site = stoi(split_str[3]);
            this->end_site = stoi(split_str[4]) - 1;
            this->n_sites = stoi(split_str[5]);
        }
        Match(int hap1, int hap2, int start_site, int end_site){
            this->hap1 = std::min(hap1, hap2);
            this->hap2 = std::max(hap1, hap2);
            this->start_site = start_site;
            this->end_site = end_site - 1;
            this->n_sites = end_site - start_site;
        }
        void display(){
            // simple display function
            if(this->len_cm > 0){
                std::cout << this->hap1 << "\t" << this->hap2 << "\t" << this->start_site << "\t" << this->end_site << "\t" << this->len_cm << std::endl;
            }
            else{
                std::cout << this->hap1 << "\t" << this->hap2 << "\t" << this->start_site << "\t" << this->end_site << std::endl;
            }
        }
        // to check if a seed is equivalent to another seed, we check to see if the haplotypes match, as well as the start and end sites
        bool operator==(const Match &other){
            
            if((this->hap1 == other.hap1) && (this->hap2 == other.hap2) && (this->start_site == other.start_site) && (this->end_site == other.end_site)){
                return true;
            }
            else{
                return false;
            }
        }
};

// apologies for the unwieldy arguments, will be fixed later
// nextStart moves backwards to find a potential extension before the beginning of a seed
// if another seed segment is directly adjacent to it moving backwards, the seed is discarded, since the seed will be included when that previous seed is extended forward
// OUTPUT: the next site if possible, else -1
int nextStart(int hap1, int hap2, int start, int max_gap, std::vector<int> &site_mapping, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, std::vector<std::vector<int>> &genotype_array);


// extendStart calls nextStart until either another seed segment is found, or you reach a point where you can no longer extend (gap greater than max_gap)
// OUTPUT: the next site if possible, else -1
int extendStart(int hap1, int hap2, int start, int max_gap, std::vector<int> &site_mapping, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, std::vector<std::vector<int>> &genotype_array);

// nextInclEnd moves forward to find a potential extension after the end of a seed
// OUTPUT: the next site for an extension (if possible)
int nextInclEnd(int hap1, int hap2, int incl_end, int max_gap, std::vector<int> &site_mapping, rateMapData &gen_map, float min_seed, float min_extend, int min_extend_markers, std::vector<std::vector<int>> &genotype_array);

// extendInclEnd calls nextInclEnd until it can no longer be extended
// OUTPUT: the last possible site after extension (if possible)
int extendInclEnd(int hap1, int hap2, int incl_end, int max_gap, std::vector<int> &site_mapping, rateMapData &gen_map, float min_seed, float min_extend, int min_extend_markers, std::vector<std::vector<int>> &genotype_array);

// simple function to match the output format of hap-IBD's original implentation
// OUTPUT: string in the format of hap-IBD
std::string hapToTskId(int hap);

// processSeed first checks to see if a seed can be extended backwards
// if not, the seed is discarded
// if so, it extends the seed backwards until it no longer can
// if the seed is extended backward, it is then extended forward (if possible)
// once the extension is done, if the segment is longer than the input threshold, the segment is written to the output file
// OUTPUT: the output string
std::string processSeed(int hap1, int hap2, int start, int incl_end, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, float min_output, std::vector<std::vector<int>> &genotype_array);

int extendBoundaryEnd(int hap1, int hap2, int incl_end, std::vector<int> &site_mapping, std::vector<std::vector<int>> &genotype_array);
int extendBoundaryStart(int hap1, int hap2, int start, std::vector<std::vector<int>> &genotype_array);


