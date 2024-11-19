#include <iostream>
#include "read_rate_map.hpp"
#include "extract_ibd_segments.hpp"
#include "vcf.hpp"
#include "RPBWT.hpp"
#include "PBWT.hpp"
#include "hap-ibd.hpp"
#include "statistics.hpp"

int main(int argc, char** argv){
    char* ts_file = argv[1];
    char* genetic_map_file = argv[2];
    float segment_threshold = std::stof(argv[3]);
    int n_threads = std::stoi(argv[4]);
    std::string minor_allele_frequency_cutoff = argv[5];
    std::string vcf_error_rate = argv[6];

    
    std::string file_string(ts_file);
    std::string stem = split(file_string, '/').back();
    std::string root = split(stem, '.')[0];

    rateMapData gen_map = readRateMap(genetic_map_file);
    // gen_map.write_interpolated_map("interpolated_map.txt");
    // gt_extraction_driver(ts_file, gen_map, segment_threshold, n_threads);
    // system("cat output* > gt_segments.txt");
    // system("rm output*");
    std::stringstream ss;
    ss << "tskit vcf --contig-id 20 " << ts_file << " > " << root << ".vcf";
    system(ss.str().c_str());
    std::string vcf_name = root + ".vcf";
    std::stringstream bcftools_command;
    std::string maf_vcf = root + "_maf" + minor_allele_frequency_cutoff + ".vcf";
    bcftools_command << "bcftools view -m2 -M2 -v snps -q " << minor_allele_frequency_cutoff << ":minor " << vcf_name << " > " << maf_vcf;
    system(bcftools_command.str().c_str());
    VCF v(maf_vcf);
    v.implant_error(std::stod(vcf_error_rate));
    std::string vcf_filtered_w_error = root + "_maf" + minor_allele_frequency_cutoff + "_e" + vcf_error_rate + ".vcf";
    v.write_to_file(vcf_filtered_w_error);
    int ret = runRPBWT(vcf_filtered_w_error, "write_to", 10000);
    ret = runPBWT(vcf_filtered_w_error, "write_to", "interpolated_map.txt", 10000, 20, 20, 1, 0.05);
    // hapIBDCpp obj("tmp.smooth.vcf", genetic_map_file, n_threads);
    // std::cout << compute_power("../../research_github/pipeline/gt_segments.txt", "output.txt") << std::endl;
    // std::cout << compute_accuracy("../../research_github/pipeline/gt_segments.txt", "output.txt") << std::endl;


    
    
}