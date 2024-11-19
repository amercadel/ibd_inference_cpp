#pragma once
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

struct interval{
    int start;
    int end;
};

void filter_hap(const std::string &input_fp, const std::string &output_fp, float threshold);
void filter_gt(const std::string &input_fp, const std::string &output_fp, float threshold);
void format_segments(const std::string &hap_ibd_fp, const std::string &gt_fp, const std::string &output_hap_ibd_fp, const std::string &output_gt_fp);
std::vector<std::string> split(const std::string &line, char delim);

// simple binary search for an insertion index if you need to only search a certain subset of a vector
// OUTPUT: where a value would be inserted into vector to keep it sorted
template<typename T>
int findInsertionIndex(std::vector<T> &vec, int from, int to, T val){
    int left = from;
    int right = to;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return left;
}

// simple binary search for an insertion index across and entire vector
// OUTPUT: where a value would be inserted into vector to keep it sorted
template<typename T>
int findInsertionIndex(std::vector<T> &vec, T val){
    int left = 0;
    int right = vec.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return left;
}


// simple binary search for finding a value, 
// OUTPUT: index if found, -1 if not
template<typename T>
int findVectorIndex(std::vector<T> &vec, T val){
    int left = 0;
    int right = vec.size() - 1;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vec[mid] == val) {
            return mid;
        } else if (vec[mid] < val) {
            left = mid + 1;
        } else {
            right = mid - 1;
        }
    }
    return -1;
}

// simple rounding function
template<typename T>
T roundToNDigits(T num, int n_digits){
    auto scale = pow(10.0, n_digits);
    return round(num * scale) / scale;
}

// generates overlapping windows to guarantee that a seed is found, but allowing for parallelization
// OUTPUT: a vector of pairs representing where to split the VCF to guarantee seed will be found
std::vector<std::pair<int, int>> overlappingWindows(std::vector<float> cm, float min_seed, int min_markers, int n_threads);
int minSites(std::vector<float> &cm_mapping, float min_seed);
