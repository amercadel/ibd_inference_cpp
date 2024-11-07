#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

struct interval{
    int start;
    int end;
};

void filter_hap(const std::string &input_fp, const std::string &output_fp, double threshold);
void filter_gt(const std::string &input_fp, const std::string &output_fp, double threshold);
void format_segments(const std::string &hap_ibd_fp, const std::string &gt_fp, const std::string &output_hap_ibd_fp, const std::string &output_gt_fp);
std::vector<std::string> split(const std::string &line, char delim);
void sort_file(const std::string &input_file, std::vector<int> sort_by, bool reverse);