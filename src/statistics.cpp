#include <statistics.hpp>
#include <segment.hpp>


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

double compute_length_accuracy(std::string gt_fp, std::string rp_fp){
    std::ifstream gt_input(gt_fp);
    std::ifstream rp_input(rp_fp);
    std::string line_g, line_r;
    std::getline(gt_input, line_g);
    std::getline(rp_input, line_r);
    IBDSegment gt_obj(line_g);
    IBDSegment rp_obj(line_r);
    std::vector<double> coverage_array;
    std::vector<IBDSegment> ibd_r;
    while(rp_obj.index1 != -1){
        while((gt_obj < rp_obj) && (gt_obj.index1 != -1)){
            std::getline(gt_input, line_g);
            gt_obj = IBDSegment(line_g);
        }
        while((gt_obj == rp_obj) && (gt_obj.index1 != -2)){
            ibd_r.push_back(gt_obj);
            std::getline(gt_input, line_g);
            gt_obj = IBDSegment(line_g);
        }
        double max_cov = 0.0;
        for(int i = 0; i < ibd_r.size(); i++){
            double cov = ibd_r[i].getCoverage(rp_obj);
            if(cov > max_cov){
                max_cov = cov;
            }
        }
        coverage_array.push_back(max_cov);
        std::getline(rp_input, line_r);
        IBDSegment rp_obj_new(line_r);
        if(!(rp_obj_new == rp_obj)){
            ibd_r.clear();
        }
        rp_obj = rp_obj_new;
    }
    gt_input.close();
    rp_input.close();
    double sum = 0;
    for(int i = 0; i < coverage_array.size(); i++){
        sum = sum + coverage_array[i];
    }
    return static_cast<double>(sum / coverage_array.size());
}

double compute_accuracy(std::string gt_fp, std::string rp_fp){
    std::ifstream gt_input(gt_fp);
    std::ifstream rp_input(rp_fp);
    std::string line_g, line_r;
    std::getline(gt_input, line_g);
    std::getline(rp_input, line_r);
    IBDSegment gt_obj(line_g);
    IBDSegment rp_obj(line_r);
    std::vector<IBDSegment> ibd_r;
    int num_covered = 0;
    int num_not_covered = 0;
    int num_v = 0;
    while(rp_obj.index1 != -1){
        while(gt_obj < rp_obj and gt_obj.index1 != -1){
            std::getline(gt_input, line_g);
            gt_obj = IBDSegment(line_g);
        }
        while (gt_obj == rp_obj and gt_obj.index1!=-1){
            ibd_r.push_back(gt_obj);
            std::getline(gt_input, line_g);
            gt_obj = IBDSegment(line_g);
        }
        double proportion_sum = 0;
        bool found = false;
        for(int i = 0; i < ibd_r.size(); i++){
            proportion_sum = ibd_r[i].getCoverage(rp_obj);
            if(proportion_sum >= 0.5){
                num_covered++;
                found = true;
                break;
            }
            if (found == false){
                num_not_covered++;
            }
            num_v++;
            std::getline(rp_input, line_r);
            IBDSegment rp_obj_new(line_r);
            if(!(rp_obj_new == rp_obj)){
                ibd_r.clear();
            }
            rp_obj = rp_obj_new;
        }
    }
    gt_input.close();
    rp_input.close();

    return static_cast<double>(num_covered / num_v);
}

double compute_power(std::string gt_fp, std::string rp_fp){
    std::ifstream gt_input(gt_fp);
    std::ifstream rp_input(rp_fp);
    std::string line_g, line_r;
    std::getline(gt_input, line_g);
    std::getline(rp_input, line_r);
    IBDSegment gt_obj(line_g);
    IBDSegment rp_obj(line_r);
    
    std::vector<IBDSegment> ibd_r;
    int num_v = 0;
    double total_sum = 0.0;
    while(gt_obj.index1 != -1){
        while(rp_obj < gt_obj && rp_obj.index1 != -1){
            std::getline(rp_input, line_r);
            rp_obj = IBDSegment(line_r);
        }
        while(rp_obj == gt_obj && rp_obj.index1 != -1){
            ibd_r.push_back(rp_obj);
            std::getline(rp_input, line_r);
            rp_obj = IBDSegment(line_r);
        }
        double proportion_sum = 0.0;
        for(int i = 0; i < ibd_r.size(); i++){
            proportion_sum = proportion_sum + ibd_r[i].getCoverage(gt_obj);
        }
        total_sum = total_sum + proportion_sum;
        num_v++;
        std::getline(gt_input, line_g);
        IBDSegment gt_obj_new(line_g);
        if(!(gt_obj_new == gt_obj)){
            ibd_r.clear();
        }
        gt_obj = gt_obj_new;
    }
    gt_input.close();
    rp_input.close();
    return static_cast<double>(total_sum / num_v);
}

double compute_accumulative_power(std::string gt_fp, std::string rp_fp){
    std::ifstream gt_input(gt_fp);
    std::ifstream rp_input(rp_fp);
    std::string line_g, line_r;
    std::getline(gt_input, line_g);
    std::getline(rp_input, line_r);
    IBDSegment gt_obj(line_g);
    IBDSegment rp_obj(line_r);
    std::vector<IBDSegment> ibd_r;
    int num_v = 0;
    double total_sum = 0.0;
    while(gt_obj.index1 != -1){
        while(rp_obj < gt_obj && rp_obj.index1 != -1){
            std::getline(rp_input, line_r);
            rp_obj = IBDSegment(line_r);
        }
        while(rp_obj == gt_obj && rp_obj.index1 != -1){
            ibd_r.push_back(rp_obj);
            std::getline(rp_input, line_r);
            rp_obj = IBDSegment(line_r);
        }
        double proportion_sum = 0.0;
        std::vector<interval> intervals;
        for(int i = 0; i < ibd_r.size(); i++){
            intervals.push_back(ibd_r[i].segment_interval);
        }
        std::vector<interval> merged_intervals = merge_intervals(intervals);
        for(int i = 0; i < merged_intervals.size(); i++){
            proportion_sum = proportion_sum + getCoverage(merged_intervals[i].start, merged_intervals[i].end, gt_obj.start, gt_obj.end);
        }
        total_sum = total_sum + proportion_sum;
        num_v++;
        std::getline(gt_input, line_g);
        IBDSegment gt_obj_new = IBDSegment(line_g);

        if(!(gt_obj_new == gt_obj)){
            ibd_r.clear();
        }
        gt_obj = gt_obj_new;
    }
    gt_input.close();
    rp_input.close();
    return static_cast<double>(total_sum / num_v);
}