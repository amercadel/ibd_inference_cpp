#pragma once
#include <utils.hpp>


std::vector<interval> merge_intervals(std::vector<interval> intervals);

float compute_length_accuracy(std::string gt_fp, std::string rp_fp);
float compute_accuracy(std::string gt_fp, std::string rp_fp);
float compute_power(std::string gt_fp, std::string rp_fp);
float compute_accumulative_power(std::string gt_fp, std::string rp_fp);

