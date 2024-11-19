
#include "tskit.h"
#include "read_rate_map.hpp"
#include "extract_ibd_segments.hpp"

std::mutex mtx;



std::vector<std::vector<int>> createMRCATable(const tsk_treeseq_t &ts){
    int n_samples = ts.num_samples;
    std::vector<std::vector<int>> mrca(n_samples, std::vector(n_samples, 0));
    return mrca;
}
std::vector<std::vector<int>> createLastLeftTable(const tsk_treeseq_t &ts){
    int n_samples = ts.num_samples;
    std::vector<std::vector<int>> last_left(n_samples, std::vector(n_samples, 0));
    return last_left;
};


void extract_segments(rateMapData &gen_map, const tsk_treeseq_t &ts, size_t start_index, size_t end_index, float min_cutoff){
    std::vector<float> &map = gen_map.interpolated_cm;
    int ret = 0;
    tsk_tree_t tree;
    check_tsk_error(ret);
    ret = tsk_tree_init(&tree, &ts, 0);
    check_tsk_error(ret);
    ret = tsk_tree_first(&tree);
    check_tsk_error(ret);
    int n_samples = ts.num_samples; // equivalent to id_end
    std::vector<std::vector<int>> mrca_last = createMRCATable(ts);
    std::vector<std::vector<int>> last_left = createLastLeftTable(ts);
    float sequence_length = tsk_treeseq_get_sequence_length(&ts);
    tsk_id_t mrca;
    std::ofstream output_file;
    char file_name[100];
    snprintf(file_name, sizeof(file_name), "output_%zu_%zu.txt", start_index, end_index);
    output_file.open(file_name, std::ofstream::out);
    if (!output_file.is_open()) {
        std::cerr << "Failed to open output file." << std::endl;
    }
    for(int i = 0; i < n_samples; i++){
        for(int j = 0; j < n_samples; j++){
            ret = tsk_tree_get_mrca(&tree, i, j, &mrca);
            check_tsk_error(ret);
            mrca_last[i][j] = mrca;
            last_left[i][j] = tree.interval.left;
        }
    }
    int tree_cnt = 0;
    int n_trees = tsk_treeseq_get_num_trees(&ts);

    int min_tree_subsample = 5000;
    int left, bp_end, bp_start, last_tree_pos;
    float gen_end, gen_start;
    for (ret = tsk_tree_first(&tree); ret == TSK_TREE_OK; ret = tsk_tree_next(&tree)){
        if (tree.interval.left - last_tree_pos < min_tree_subsample){
            tree_cnt++;
            if ((tree_cnt % 50000 == 0)){
            std::cout << tree_cnt << " out of " << n_trees << " visited" << std::endl;
            }

            continue;
        }
        last_tree_pos = tree.interval.left;
        for(size_t i = start_index; i < end_index; i++){
            for(size_t j = i + 1; j < n_samples; j++){
                ret = tsk_tree_get_mrca(&tree, i, j, &mrca);
                if(mrca_last[i][j] != mrca){
                    left = tree.interval.left;
                    bp_end = left;
                    bp_start = last_left[i][j];
                    if (bp_end > map.size()){
                        gen_end = map[map.size() - 1];
                    }else{
                        gen_end = map[bp_end];
                    }
                    if (bp_start > map.size()){
                        gen_start = map[map.size()];
                    }else{
                        gen_start = map[bp_start];
                    }
                    if (gen_end - gen_start < min_cutoff){
                        last_left[i][j] = left;
                        mrca_last[i][j] = mrca;
                        continue;
                    }
                    if (gen_end - gen_start > min_cutoff){
                        float len = gen_end - gen_start;
                        ibd_segment seg = ibd_segment(i, j, bp_start, bp_end, gen_start, gen_end, len);
                        std::string out = seg.to_string();
                        output_file << out << std::endl;
                        last_left[i][j] = left;
                        mrca_last[i][j] = mrca;
                    }

                }

            }
        }
        tree_cnt++;
        if ((tree_cnt % 50000 == 0)){
            std::cout << tree_cnt << " out of " << n_trees << " visited" << std::endl;
        }


    }

    for (size_t i = start_index; i < end_index; i++){
        for (size_t j = i + 1; j < n_samples; j++){
            bp_start = last_left[i][j];
            
            bp_end = static_cast<int>(sequence_length);
            gen_end = 0;
            gen_start = 0;
            if (bp_end > map.size()){
                gen_end = map[map.size() - 1];
            }else{
                gen_end = map[bp_end];   
            }
            if(bp_start > map.size()){
                gen_start = map[map.size() - 1];

            }else{
                gen_start = map[bp_start];
            }
            
            if (gen_end - gen_start > min_cutoff){
                float len = gen_end - gen_start;
                ibd_segment seg = ibd_segment(i, j, bp_start, bp_end, gen_start, gen_end, len);
                std::string out = seg.to_string();
                output_file << out << std::endl;
            }
        }
    }
    output_file.close();

    tsk_tree_free(&tree);
}

std::vector<std::pair<int, int>> generate_cutoffs(int n_samples, int n_threads){
    int total_comparisons = n_samples * (n_samples - 1) / 2;
    int comparisons_per_thread = total_comparisons / n_threads;
    std::vector<std::pair<int, int>> cutoffs;
    size_t start = 0;
    int accumulated_comparisons = 0;
    for (size_t i = 0; i < n_threads; ++i) {
        size_t end = start;

        int thread_comparisons = 0;
        while (end < n_samples - 1 && thread_comparisons < comparisons_per_thread) {
            int comparisons_for_index = n_samples - 1 - end;
            if (thread_comparisons + comparisons_for_index > comparisons_per_thread) {
                break;
            }
            thread_comparisons += comparisons_for_index;
            end++;
        }
        std::pair<int, int> p;
        p.first = start;
        if(i == n_threads - 1){
            p.second = n_samples;
        }
        else{
            p.second = end;
        }
        accumulated_comparisons += thread_comparisons;
        start = end;
        cutoffs.push_back(p);
    }
    return cutoffs;
}

void gt_extraction_driver(char* ts_file, rateMapData gen_map, float minimum_cutoff, int n_threads){
    tsk_treeseq_t ts;
    int ret = 0;
    ret = tsk_treeseq_load(&ts, ts_file, 0);
    check_tsk_error(ret);
    std::thread threads[n_threads];
    std::vector<std::pair<int, int>> cutoffs = generate_cutoffs(ts.num_samples, n_threads);
    for(size_t i = 0; i < n_threads; i++){
        threads[i] = std::thread(extract_segments, std::ref(gen_map), std::ref(ts), cutoffs[i].first, cutoffs[i].second, minimum_cutoff);
    }
    for (auto& th : threads) th.join();
    tsk_treeseq_free(&ts);
}