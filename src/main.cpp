#include <iostream>
#include <chrono>
#include "read_rate_map.hpp"
#include "extract_ibd_segments.hpp"
#include "vcf.hpp"
#include "RPBWT.hpp"
#include "PBWT.hpp"
#include "hap-ibd.hpp"
#include "statistics.hpp"

int main(int argc, char** argv){
    auto start_main = std::chrono::high_resolution_clock::now();
    char* ts_file = argv[1];
    char* genetic_map_file = argv[2];
    float segment_threshold = std::stof(argv[3]);
    int n_threads = std::stoi(argv[4]);
    std::string minor_allele_frequency_cutoff = argv[5];
    std::string vcf_error_rate = argv[6];

    
    std::string file_string(ts_file);
    std::string stem = split(file_string, '/').back();
    std::string root = split(stem, '.')[0];

    auto start = std::chrono::high_resolution_clock::now();
    rateMapData gen_map = readRateMap(genetic_map_file);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time taken to read rate map: " << elapsed.count() << " seconds" << std::endl;

    auto start_interpolation = std::chrono::high_resolution_clock::now();
    gen_map.write_interpolated_map("interpolated_map.txt");
    auto end_interpolation = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_interpolation = end_interpolation - start_interpolation;
    std::cout << "Time taken to write interpolated map: " << elapsed_interpolation.count() << " seconds" << std::endl;


    auto start_extraction = std::chrono::high_resolution_clock::now();
    gt_extraction_driver(ts_file, gen_map, segment_threshold, n_threads);
    auto end_extraction = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_extraction = end_extraction - start_extraction;
    std::cout << "Time taken to extract GT segments: " << elapsed_extraction.count() << " seconds" << std::endl;


    system("cat output* > gt_segments.txt");
    system("rm output*");

    std::stringstream ss;
    ss << "tskit vcf --contig-id 20 " << ts_file << " > " << root << ".vcf";
    auto start_vcf = std::chrono::high_resolution_clock::now();
    system(ss.str().c_str());
    auto end_vcf = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_vcf = end_vcf - start_vcf;
    std::cout << "Time taken to generate VCF: " << elapsed_vcf.count() << " seconds" << std::endl;

    std::string vcf_name = root + ".vcf";
    std::stringstream bcftools_command;
    std::string maf_vcf = root + "_maf" + minor_allele_frequency_cutoff + ".vcf";

    bcftools_command << "bcftools view -m2 -M2 -v snps -q " << minor_allele_frequency_cutoff << ":minor " << vcf_name << " > " << maf_vcf;
    auto start_bcftools = std::chrono::high_resolution_clock::now();
    system(bcftools_command.str().c_str());
    auto end_bcftools = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_bcftools = end_bcftools - start_bcftools;
    std::cout << "Time taken to filter VCF with bcftools: " << elapsed_bcftools.count() << " seconds" << std::endl;

    auto start_vcf_read = std::chrono::high_resolution_clock::now();
    VCF v(maf_vcf);
    auto end_vcf_read = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_vcf_read = end_vcf_read - start_vcf_read;
    std::cout << "Time taken to read VCF: " << elapsed_vcf_read.count() << " seconds" << std::endl;

    auto start_implant_error = std::chrono::high_resolution_clock::now();
    v.implant_error(std::stod(vcf_error_rate));
    auto end_implant_error = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_implant_error = end_implant_error - start_implant_error;
    std::cout << "Time taken to implant error: " << elapsed_implant_error.count() << " seconds" << std::endl;

    std::string vcf_filtered_w_error = root + "_maf" + minor_allele_frequency_cutoff + "_e" + vcf_error_rate + ".vcf";
    auto start_write_to_file = std::chrono::high_resolution_clock::now();
    v.write_to_file(vcf_filtered_w_error);
    auto end_write_to_file = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_write_to_file = end_write_to_file - start_write_to_file;
    std::cout << "Time taken to write VCF to file: " << elapsed_write_to_file.count() << " seconds" << std::endl;

    auto start_pbwt = std::chrono::high_resolution_clock::now();
    int ret = runRPBWT(maf_vcf, "write_to", 10000);
    ret = runPBWT(vcf_filtered_w_error, "write_to", "interpolated_map.txt", 10000, 20, 20, 1, 0.05);
    auto end_pbwt = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_pbwt = end_pbwt - start_pbwt;
    std::cout << "Time taken to run RPBWT and PBWT: " << elapsed_pbwt.count() << " seconds" << std::endl;


    auto start_hapibd = std::chrono::high_resolution_clock::now();
    hapIBDCpp obj("write_to.smooth.vcf", genetic_map_file, n_threads);
    auto end_hapibd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_hapibd = end_hapibd - start_hapibd;
    std::cout << "Time taken to initialize hapIBDCpp: " << elapsed_hapibd.count() << " seconds" << std::endl;
    
    
    auto start_all_computations = std::chrono::high_resolution_clock::now();

    float s = compute_power("gt_segments.txt", "output.txt");
    s = compute_accuracy("gt_segments.txt", "output.txt");
    s = compute_length_accuracy("gt_segments.txt", "output.txt");
    s = compute_accumulative_power("gt_segments.txt", "output.txt");

    auto end_all_computations = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_all_computations = end_all_computations - start_all_computations;
    std::cout << "Time taken to run all computations: " << elapsed_all_computations.count() << " seconds" << std::endl;
    
    auto end_main = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_main = end_main - start_main;
    std::cout << "Total time taken by main function: " << elapsed_main.count() << " seconds" << std::endl;


    
    
}