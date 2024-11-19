#include <iostream>
#include <fstream>
#include <sstream>

#include "utils.hpp"
#include "read_rate_map.hpp"





std::vector<float> rateMapData::interpolateVector(std::vector<int> &sites, std::vector<int> &bp, std::vector<float> &cm) {
    std::vector<float> interpolated_cm;
    int n = sites.size();
    int m = bp.size();
    int j = 0;
    for (int i = 0; i < n; i++){
        while(j < m && bp[j] < sites[i]){
            j++;
        }
        if (j == 0){
            interpolated_cm.push_back(cm[0]);
        }else if (j == m){
            interpolated_cm.push_back(cm[m - 1]);
        }else{
            int x0 = bp[j - 1];
            int x1 = bp[j];
            float y0 = cm[j - 1];
            float y1 = cm[j];
            int x = sites[i];
            float y = ((y0 * (x1 - x)) + (y1 * (x - x0))) / (x1 - x0);
            interpolated_cm.push_back(y);
        }
    }
    return interpolated_cm;
}


rateMapData readRateMap(char* filename){
    std::ifstream inputFile;
    std::string line;
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    inputFile.open(filename);
    std::cout << "opening file\n";
    rateMapData res;
    while (getline(inputFile, line)) {
        int bp;
        float cm;
        std::string tmp;
        std::stringstream input_str(line); 
        std::vector<std::string> data = split(line, ' ');
        bp = std::stoi(data[3]);
        cm = std::stod(data[2]);
        bp_vec.push_back(bp);
        cm_vec.push_back(cm);

    }
    int min_tree_subsample = INT_MAX;
    for (int a = 1; a < bp_vec.size() - 1; a++){
        int sub_samp = bp_vec[a] - bp_vec[a - 1];
        if(sub_samp < min_tree_subsample){
            min_tree_subsample = sub_samp;
        }
    }
    res.bp_vec = bp_vec;
    res.cm_vec = cm_vec;
    res.min_tree_subsample = min_tree_subsample;
    res.last_bp = bp_vec[bp_vec.size() - 1];
    res.last_cm = cm_vec[cm_vec.size() - 1];
    for (int i = 0; i <= res.last_bp + 1; i ++){
        res.all_sites.push_back(i);
    }
    res.interpolated_cm = res.interpolateVector(res.all_sites, res.bp_vec, res.cm_vec);
    return res;
}

void rateMapData::write_interpolated_map(std::string output_file){
    std::ofstream out(output_file);
    for(int i = 0; i < this->interpolated_cm.size(); i++){
        std::stringstream oss;
        oss << i << "\t" << interpolated_cm[i] << "\n";
        out << oss.str();
    }
    out.close();
}


float getGeneticPosition(std::vector<float> &interpolated_cm, int site_idx){
    // assert(site_idx < interpolated_cm.size());
    // if(site_idx < interpolated_cm.size()){
    //     return interpolated_cm[site_idx];
    // }
    // else{
    //     std::cerr << "Site outside of interolated map\n";
    //     exit(EXIT_FAILURE);
    // }
    return interpolated_cm[site_idx];
	
}
std::vector<float> rateMapData::interpolateVector(std::vector<int> &sites) {
    std::vector<float> interpolated_cm;
    for(size_t c = 0; c < sites.size(); c++){
        interpolated_cm.push_back(genPos(sites[c]));
    }
    return interpolated_cm;
}

float rateMapData::genPos(int site){
    int map_size_m1 = bp_vec.size() - 1;
    float min_end_cm_dist = 5.0f;
    auto it = std::lower_bound(bp_vec.begin(), bp_vec.end(), site);
    int index = it - bp_vec.begin();
    if (it != bp_vec.end() && *it == site){
        return cm_vec[index];
    }else{
        int insertion_point = index;
        int a_index = insertion_point - 1;
        int b_index = insertion_point;
        if (a_index == map_size_m1){
            auto gp_it = std::lower_bound(cm_vec.begin(), cm_vec.end(), cm_vec.back() - min_end_cm_dist);
            insertion_point = gp_it - cm_vec.begin();
            if(insertion_point < 0){
                insertion_point = -insertion_point - 2;
            }
            a_index = std::max(insertion_point, 0);
            b_index = map_size_m1;
        }
        else if(b_index == 0){
            auto gp_it = std::lower_bound(cm_vec.begin(), cm_vec.end(), cm_vec.front() + min_end_cm_dist);
            insertion_point = gp_it - cm_vec.begin();
            if(insertion_point < 0){
                insertion_point = -insertion_point - 2;
            }
            a_index = 0;
            b_index = std::min(insertion_point, map_size_m1);
        }
        int x = site;
        int a = bp_vec[a_index];
        int b = bp_vec[b_index];
        float fa = cm_vec[a_index];
        float fb = cm_vec[b_index];
        return fa + (((float)(x - a) / (float)(b - a)) * (fb - fa));
    }
    
}

float rateMapData::interpolateBasePairToGenPos(int site){
    int idx = findVectorIndex(bp_vec, site);
    if (idx > 0){
        return cm_vec[idx];
    }
    else{
        idx = findInsertionIndex(bp_vec, site);
        if (idx == 0){
            return cm_vec[0];
        }
        else if(idx >= bp_vec.size()){
            float slope = ((cm_vec.back() - cm_vec.front()) / (bp_vec.back() - bp_vec.front()));
            float y_int = cm_vec.back() - (slope * bp_vec.back());
            return slope * site + y_int;
        }

        else{
            int x0 = bp_vec[idx - 1];
            int x1 = bp_vec[idx];
            float y0 = cm_vec[idx - 1];
            float y1 = cm_vec[idx];
            float y = ((y0 * (x1 - site)) + (y1 * (site - x0))) / (x1 - x0);
            return y;

        }
    }
}


rateMapData readRateMap(char* filename, std::vector<int> &sites){
    std::ifstream inputFile;
    std::string line;
    std::vector<int> bp_vec;
    std::vector<float> cm_vec;
    inputFile.open(filename);
    std::cout << "reading genetic map file\n";
    rateMapData res;
    while (getline(inputFile, line)) {
        int bp;
        float cm;
        std::string tmp;
        std::stringstream input_str(line); 
        std::vector<std::string> data = split(line, ' ');
        bp = std::stoi(data[3]);
        cm = std::stod(data[2]);
        bp_vec.push_back(bp);
        cm_vec.push_back(cm);

    }
    std::cout << "rate map loaded into memory\n";
    res.bp_vec = bp_vec;
    res.cm_vec = cm_vec;
    res.last_bp = bp_vec[bp_vec.size() - 1];
    res.last_cm = cm_vec[cm_vec.size() - 1];
    std::cout << "interpolating rate map\n";
    res.interpolated_cm = res.interpolateVector(sites);
    std::cout << "rate map interpolated\n";
    return res;
}

