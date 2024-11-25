#include <statistics.hpp>
#include <segment.hpp>
#include <iostream>


std::vector<interval> merge_intervals(std::vector<interval> intervals){
    if(intervals.size() <= 1){
        return intervals;
    }
    std::vector<interval> merged;
    for(int i = 0; i < intervals.size(); i++){
        if((merged.size() == 0) || merged.back().end < intervals[i].start){
            merged.push_back(intervals[i]);
        }
        else{
            merged.back().end = std::max(merged.back().end, intervals[i].end);
        }
    }
    return merged;
}

float compute_length_accuracy(std::string gt_fp, std::string rp_fp){
    std::ifstream gt_input(gt_fp);
    std::ifstream rp_input(rp_fp);
    std::string line_g, line_r;
    std::vector<IBDSegment> ground_truth;
    std::vector<IBDSegment> reported;
    std::vector<float> cov_array;
    while (std::getline(gt_input, line_g)) {
        ground_truth.push_back(IBDSegment(line_g));
    }
    while (std::getline(rp_input, line_r)) {
        reported.push_back(IBDSegment(line_r));
    }
    for(int i = 0; i < reported.size(); i++){
        float max_cov = 0;
        for(int j = 0; j < ground_truth.size(); j++){
            if (!(reported[i] == ground_truth[j])){
                continue;
            }
            float cov = reported[i].getCoverage(ground_truth[j]);
            if (cov > max_cov){
                max_cov = cov;
            }
        }
        cov_array.push_back(max_cov);
    }
    float sum_cov = 0;
    for (float cov : cov_array) {
        sum_cov += cov;
    }
    return sum_cov / cov_array.size();
    
}

float compute_accuracy(std::string gt_fp, std::string rp_fp){
    std::ifstream gt_input(gt_fp);
    std::ifstream rp_input(rp_fp);
    std::string line_g, line_r;
    std::vector<IBDSegment> ground_truth;
    std::vector<IBDSegment> reported;
    std::vector<float> cov_array;
    while (std::getline(gt_input, line_g)) {
        ground_truth.push_back(IBDSegment(line_g));
    }
    while (std::getline(rp_input, line_r)) {
        reported.push_back(IBDSegment(line_r));
    }
    int n_covered = 0;
    
    for(int i = 0; i < reported.size(); i++){
        for(int j = 0; j < ground_truth.size(); j++){
            if(reported[i].index1 == ground_truth[j].index1 && reported[i].index2 == ground_truth[j].index2){
                float cov = reported[i].getCoverage(ground_truth[j]);
                if(cov > 0.5f){
                    n_covered++;
                    break;
                }
            }
        }
    }
    return static_cast<float>(n_covered) / reported.size();

}

float compute_power(std::string gt_fp, std::string rp_fp){
    std::ifstream gt_input(gt_fp);
    std::ifstream rp_input(rp_fp);
    std::string line_g, line_r;
    std::vector<IBDSegment> ground_truth;
    std::vector<IBDSegment> reported;
    std::vector<float> cov_array;
    while (std::getline(gt_input, line_g)) {
        ground_truth.push_back(IBDSegment(line_g));
    }
    while (std::getline(rp_input, line_r)) {
        reported.push_back(IBDSegment(line_r));
    }
    
    for(int i = 0; i < ground_truth.size(); i++){
        float max_cov = 0.0f;
        for(int j = 0; j < reported.size(); j++){
            if (ground_truth[i].index1 == reported[j].index1 && ground_truth[i].index2 == reported[j].index2){
                float cov = ground_truth[i].getCoverage(reported[j]);
                if (cov > max_cov){
                    max_cov = cov;
                }
            }
        }
        cov_array.push_back(max_cov);
    }
    float sum_cov = 0;
    for (float cov : cov_array) {
        sum_cov += cov;
    }
    return static_cast<float>(sum_cov / cov_array.size());

    
}

float compute_accumulative_power(std::string gt_fp, std::string rp_fp){
    std::ifstream gt_input(gt_fp);
    std::ifstream rp_input(rp_fp);
    std::string line_g, line_r;
    std::vector<IBDSegment> ground_truth;
    std::vector<IBDSegment> reported;
    std::vector<float> cov_array;
    while (std::getline(gt_input, line_g)) {
        ground_truth.push_back(IBDSegment(line_g));
    }
    while (std::getline(rp_input, line_r)) {
        reported.push_back(IBDSegment(line_r));
    }
    for(int i = 0; i < ground_truth.size(); i++){
        float proportion = 0.0f;
        for(int j = 0; j < reported.size(); j++){
            if (!(ground_truth[i] == reported[j])){
                continue;
            }
            float cov = ground_truth[i].getCoverage(reported[j]);
            proportion = proportion + cov;
        }
        cov_array.push_back(proportion);
    }
    float sum_cov = 0;
    for (float cov : cov_array) {
        sum_cov += cov;
    }
    return sum_cov / cov_array.size();

}