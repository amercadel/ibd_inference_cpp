#pragma once
#include <vector>
#include <iostream>
#include <string>
#include <thread>
#include "tskit.h"

#define check_tsk_error(val)                                                            \
    if (val < 0) {                                                                      \
        fprintf(stderr, "line %d: %s", __LINE__, tsk_strerror(val));                    \
        exit(EXIT_FAILURE);                                                             \
    }


struct ibd_segment{
    int id1;
    int id2;
    int bp_start;
    int bp_end;
    float cm_start;
    float cm_end;
    float len_cm;
    void display(){
        std::cout << id1 << "\t" << id2 << "\t" << bp_start << "\t" << bp_end << "\t" << cm_start << "\t" << \
        cm_end << "\t" << len_cm << std::endl;
        }
    
    std::string to_string() {
        return std::to_string(id1) + "\t" + std::to_string(id2) + "\t" + std::to_string(bp_start) + "\t" +
                std::to_string(bp_end) + "\t" + std::to_string(cm_start) + "\t" + std::to_string(cm_end) + "\t" +
                std::to_string(len_cm);
    }
    
    ibd_segment(int id1, int id2, int bp_start, int bp_end, float cm_start, float cm_end, float len_cm) :
        id1(id1), id2(id2), bp_start(bp_start), bp_end(bp_end), cm_start(cm_start), cm_end(cm_end), len_cm(len_cm) {}

};

std::vector<std::vector<int>> createMRCATable(const tsk_treeseq_t &ts);
std::vector<std::vector<int>> createLastLeftTable(const tsk_treeseq_t &ts);
void extract_segments(rateMapData &gen_map, const tsk_treeseq_t &ts, size_t start_index, size_t end_index, float min_cutoff);
std::vector<std::pair<int, int>> generate_cutoffs(int n_samples, int n_cores);
void gt_extraction_driver(char* ts_file, rateMapData gen_map, float minimum_cutoff, int n_threads);