#pragma once
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <unordered_map>
#include "utils.hpp"


class VCF{
    public:
        std::string file_path;
        std::vector<std::string> header_info;
        std::vector<std::vector<std::string>> all_data;
        int n_sites;
        int n_individuals;
        int n_samples;
        std::vector<int> sites;
        std::unordered_map<std::string, int> ids_map;
        std::vector<std::string> individuals;
        VCF(const std::string& file_path){
            this->file_path = file_path;
            std::ifstream fp(this->file_path);
            std::string line;
            while(std::getline(fp, line)){
                if (line.at(0) == '#'){
                    header_info.push_back(line);
                }
                else{
                    std::vector<std::string> tmp = split(line, '\t');
                    this->sites.push_back(std::stoi(tmp[1]));
                    all_data.push_back(tmp);
                }
            }
            this->n_sites = this->all_data.size();
            this->n_individuals = this->all_data[0].size() - 9;
            this->n_samples = this->n_individuals * 2;
            std::vector<std::string> tmp = split(this->header_info.back(), '\t');
            for(int i = 9; i < tmp.size(); i++){
                this->ids_map[tmp[i]] = i + 9;
            }
        }
        void implant_error(double error_rate);

        void write_to_file(std::string output_file);




    private:
        std::string flip_allele(const std::string allele);

};

int count_p_smoother_corrections(VCF unsmooth_vcf, VCF smooth_vcf);