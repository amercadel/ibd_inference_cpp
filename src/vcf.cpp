#include "vcf.hpp"


std::string VCF::flip_allele(const std::string allele){
    if(allele == "0"){
        return "1";
    }
    else{
        return "0";
    }
}

void VCF::implant_error(double error_rate){
    std::random_device rd;
    for(int i = 9; i < this->all_data[0].size(); i++){
        
        int n_sites_flipped = this->n_sites * error_rate;
        std::vector<int> range(this->n_sites);
        for (int c = 0; c < this->n_sites; c++) {
            range[c] = c;
        }
        std::mt19937 g(rd());
        std::shuffle(range.begin(), range.end(), g);
        range.resize(n_sites_flipped);
        for(int j = 0; j < range.size(); j++){
            int allele_to_flip = rand() % 2;
            std::vector<std::string> alleles = split(this->all_data[range[j]][i], '|');
            if(allele_to_flip == 0){
                alleles[0] = flip_allele(alleles[0]);
            }
            else{
                alleles[1] = flip_allele(alleles[1]);
            }
            std::stringstream tmp;
            tmp << alleles[0] << "|" << alleles[1];
            all_data[range[j]][i] = tmp.str();
        }
    }
}
void VCF::write_to_file(std::string output_file){
    std::ofstream out(output_file);
    for(int i = 0; i < this->header_info.size(); i++){
        out << this->header_info[i] << "\n";
    }
    for(int i = 0; i < this->all_data.size(); i++){
        std::stringstream tmp;
        for(int j = 0; j < this->all_data[i].size() - 1; j++){
            tmp << all_data[i][j] << "\t";
        }
        tmp << all_data[i].back() << "\n";
        out << tmp.str();
    }
    out.close();
}

int count_p_smoother_corrections(VCF unsmooth_vcf, VCF smooth_vcf){
    int n_corrections = 0;
    for(int i = 0; i < unsmooth_vcf.all_data.size(); i++){
        for(int j = 9; j < unsmooth_vcf.all_data[0].size(); j++){
            std::vector<std::string> unsmooth_vals = split(unsmooth_vcf.all_data[i][j], '|');
            std::vector<std::string> smooth_vals = split(smooth_vcf.all_data[i][j], '|');
            if(unsmooth_vals[0] != smooth_vals[0]){
                n_corrections++;
            }
            if(unsmooth_vals[1] != smooth_vals[1]){
                n_corrections++;
            }
        }
    }
    return n_corrections;
}
