#include <cassert>
#include "match.hpp"
#include "utils.hpp"

int nextStart(int hap1, int hap2, int start, int max_gap, std::vector<int> &site_mapping, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, std::vector<std::vector<int>> &genotype_array){
    if (start < 2 || max_gap < 0){
        return start;
    }
    int m = start - 1;
    int first_mismatch_pos = site_mapping[m];

    int first_match = start - 2;
    while(m > 0){
        --m;
        int a1 = genotype_array[m][hap1];
        int a2 = genotype_array[m][hap2];
        // std::cout << m << " " << a1 << " " << a2 << std::endl;
        if(a1!=a2){
            
            if((first_mismatch_pos - site_mapping[m]) > max_gap){
                ++m;
                break;
            }
            else if (m > 0){
                first_match = m - 1;
            }
        }
    }
    
    float len = gen_map.interpolated_cm[first_match] - gen_map.interpolated_cm[m];
    if (len >= min_seed && ((first_match - m) >= min_seed_markers - 1)){
        return -1;
    }
    else{
        int ret = (len < min_extend || ((first_match - m) < min_extend_markers - 1)) ? start : m;
        return ret;
    }
}


int extendStart(int hap1, int hap2, int start, int max_gap, std::vector<int> &site_mapping, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, std::vector<std::vector<int>> &genotype_array){
    int prev_start = start;
    int next_start = nextStart(hap1, hap2, prev_start, max_gap, site_mapping, gen_map, min_seed, min_extend, min_seed_markers, min_extend_markers, genotype_array);
    while (next_start>=0 && (next_start < prev_start)) {
        prev_start = next_start;
        next_start = nextStart(hap1, hap2, prev_start, max_gap, site_mapping, gen_map, min_seed, min_extend, min_seed_markers, min_extend_markers, genotype_array);
    }
    return next_start;
}

int nextInclEnd(int hap1, int hap2, int incl_end, int max_gap, std::vector<int> &site_mapping, rateMapData &gen_map, float min_seed, float min_extend, int min_extend_markers, std::vector<std::vector<int>> &genotype_array){
    int last_marker = site_mapping.size() - 1;
    if ((incl_end>(last_marker - 2)) || (max_gap < 0)) {
        return incl_end;
    }
    int m = incl_end + 1;
    int first_mismatch_pos = site_mapping[m];
    int first_match = incl_end + 2;
    while(m < last_marker){
        ++m;
        
        int a1 = genotype_array[m][hap1];
        int a2 = genotype_array[m][hap2];
        if(a1 != a2){
            if ((site_mapping[m] - first_mismatch_pos) > max_gap){
                --m;
                break;
            }
            else if(m < last_marker){
                first_match = m + 1;
            }
        }

    }
    float len = (gen_map.interpolated_cm[m] - gen_map.interpolated_cm[first_match]);
    return (len<min_extend || (m-first_match)<(min_extend_markers - 1)) ? incl_end : m;
}

int extendInclEnd(int hap1, int hap2, int incl_end, int max_gap, std::vector<int> &site_mapping, rateMapData &gen_map, float min_seed, float min_extend, int min_extend_markers, std::vector<std::vector<int>> &genotype_array){
    int last_marker = site_mapping.size() - 1;
    while (incl_end<last_marker && (genotype_array[incl_end + 1][hap1] == genotype_array[incl_end + 1][hap2])) {
        ++incl_end;
    }
    int prev_incl_end = incl_end;
    int next_incl_end = nextInclEnd(hap1, hap2, prev_incl_end, max_gap, site_mapping, gen_map, min_seed, min_extend, min_extend_markers, genotype_array);
    while (next_incl_end>prev_incl_end) {
        prev_incl_end = next_incl_end;
        next_incl_end = nextInclEnd(hap1, hap2, prev_incl_end, max_gap, site_mapping, gen_map, min_seed, min_extend, min_extend_markers, genotype_array);
    }
    return next_incl_end;
}


std::string hapToTskId(int hap){
    int id = hap / 2;
    int hap_id = (hap % 2) + 1;
    std::string ret;
    ret = "tsk_" + std::to_string(id) + "\t" + std::to_string(hap_id);
    return ret;
}


std::string processSeed(int hap1, int hap2, int start, int incl_end, int max_gap, std::vector<int> &site_mapping, std::vector<Match> &matches, rateMapData &gen_map, float min_seed, float min_extend, int min_seed_markers, int min_extend_markers, float min_output, std::vector<std::vector<int>> &genotype_array){
    std::stringstream out;
    
    start = extendStart(hap1, hap2, start, max_gap, site_mapping, gen_map, min_seed, min_extend, min_seed_markers, min_extend_markers, genotype_array);
    if (start>=0) {
        incl_end = extendInclEnd(hap1, hap2, incl_end, max_gap, site_mapping, gen_map, min_seed, min_extend, min_extend_markers, genotype_array);
        if (incl_end >= site_mapping.size()){
            incl_end = site_mapping.size() - 1;
        }
        if ((gen_map.interpolated_cm[incl_end] - gen_map.interpolated_cm[start])>=min_output) {
            out << hapToTskId(hap1) << "\t" << hapToTskId(hap2) << "\t" << "20" << "\t" << site_mapping[start] << "\t" << site_mapping[incl_end] << "\t" << roundToNDigits(gen_map.interpolated_cm[incl_end] - gen_map.interpolated_cm[start], 3) << "\n";
            }
        return out.str();
        }
    else{
        return "";
    }
}

int extendBoundaryStart(int hap1, int hap2, int start, std::vector<std::vector<int>> &genotype_array){
    if (start == 0){
        return start;
    }
    int m = start - 1;
    int a1 = genotype_array[m][hap1];
    int a2 = genotype_array[m][hap2];
    if (a1 != a2){
        return start;
    }
    while(((a1 == a2) && (m >= 0))){
        --m;
        a1 = genotype_array[m][hap1];
        a2 = genotype_array[m][hap2];
        if(m == 0){
            return 0;
        }
    }
    return m + 1;
    
}

int extendBoundaryEnd(int hap1, int hap2, int incl_end, std::vector<int> &site_mapping, std::vector<std::vector<int>> &genotype_array){
    int last_marker = site_mapping.size() - 1;
    if (incl_end == last_marker){
        return incl_end;
    }
    int m = incl_end + 1;
    int a1 = genotype_array[m][hap1];
    int a2 = genotype_array[m][hap2];
    if (a1 != a2){
        return incl_end;
    }
    while(((a1 == a2) && (m < last_marker))){
        ++m;
        a1 = genotype_array[m][hap1];
        a2 = genotype_array[m][hap2];
    }
    return m - 1;
}







