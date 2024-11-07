#pragma once
#include <utils.hpp>


std::vector<interval> merge_intervals(std::vector<interval> intervals);

double compute_length_accuracy(std::string gt_fp, std::string rp_fp);
double compute_accuracy(std::string gt_fp, std::string rp_fp);
double compute_power(std::string gt_fp, std::string rp_fp);
double compute_accumulative_power(std::string gt_fp, std::string rp_fp);

