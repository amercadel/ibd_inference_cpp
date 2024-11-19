#include "utils.hpp"
#include "segment.hpp"


void filter_hap(const std::string &input_fp, const std::string &output_fp, float threshold){
    std::ifstream input_file(input_fp);
    std::ofstream output_file(output_fp);
    std::string line;
    std::vector<std::string> tokens;
    while(std::getline(input_file, line)){
        tokens = split(line, '\t');
        if (std::stod(tokens[7]) >= threshold){
            output_file << line;
        }
    }
    input_file.close();
    output_file.close();

}

void filter_gt(const std::string &input_fp, const std::string &output_fp, float threshold){
    std::ifstream input_file(input_fp);
    std::ofstream output_file(output_fp);
    std::string line;
    std::vector<std::string> tokens;
    while(std::getline(input_file, line)){
        tokens = split(line, '\t');
        if (std::stod(tokens[4]) >= threshold){
            output_file << line;
        }
    }
    input_file.close();
    output_file.close();
}

void format_segments(const std::string &hap_ibd_fp, const std::string &gt_fp, const std::string &output_hap_ibd_fp, const std::string &output_gt_fp){
    std::ifstream hap_ibd_input_file(hap_ibd_fp);
    std::ifstream gt_input_file(gt_fp);
    std::ofstream hap_ibd_output_file(output_hap_ibd_fp);
    std::ofstream gt_output_file(output_gt_fp);
    std::string line;
    std::vector<std::string> tokens;
    while(std::getline(hap_ibd_input_file, line)){
        tokens = split(line, '\t');
        int id1 = std::stoi(split(tokens[0], '_')[1]);
        int hap1 = std::stoi(tokens[1]) - 1;
        int id2 = std::stoi(split(tokens[2], '_')[1]);
        int hap2 = std::stoi(tokens[3]) - 1;
        std::string s = tokens[5];
        std::string e = tokens[6];
        std::string len = tokens[7];
        std::stringstream oss;
        oss << (id1 * 2 - hap1) << "\t" << (id2 * 2 - hap2) << "\t" << s << "\t" << e << "\t" << len << "\n";
        hap_ibd_output_file << oss.str();
    }
    while(std::getline(gt_input_file, line)){
        gt_output_file << line;
        
    }
    hap_ibd_input_file.close();
    hap_ibd_output_file.close();
    gt_input_file.close();
    gt_output_file.close();
}



std::vector<std::string> split(const std::string &line, char delim){
    std::vector<std::string> split_string;
    std::string str;
    std::stringstream ss(line);
    while (std::getline(ss, str, delim)){
        split_string.push_back(str);
    }
    return split_string;

}

void filter_by_overlap(const std::string &hap_ibd_fp, const std::string &gt_fp, float cutoff){
    std::ifstream hap_input(hap_ibd_fp);
    std::ifstream gt_input(gt_fp);
    std::stringstream false_ibd_fp;
    false_ibd_fp << "false_ibds_" << cutoff << ".txt";
    std::stringstream true_ibd_fp;
    true_ibd_fp << "reported_ibds_" << cutoff << ".txt";
    std::ofstream false_output(false_ibd_fp.str());
    std::ofstream true_output(true_ibd_fp.str());
    std::vector<float> coverage_array;
    std::vector<IBDSegment> ibd_r;
    std::string line_g, line_r;
    std::getline(gt_input, line_g);
    std::getline(hap_input, line_r);
    IBDSegment gt_obj(line_g);
    IBDSegment rp_obj(line_r);
    
    while(rp_obj.index1 != -1){
        while((gt_obj < rp_obj) && (gt_obj.index1 != -1)){
            if(!std::getline(gt_input, line_g)) break;
            gt_obj = IBDSegment(line_g);
        }
        while((gt_obj == rp_obj) && (gt_obj.index1 != -1)){
            ibd_r.push_back(gt_obj);
            if(!std::getline(gt_input, line_g)) break;
            gt_obj = IBDSegment(line_g);
        }
        float max_cov = 0;
        for(auto &ir: ibd_r){
            float cov = ir.getCoverage(rp_obj);
            if(cov > max_cov){
                max_cov = cov;
            }
        }
        coverage_array.push_back(max_cov);
        if(max_cov >= cutoff){
            true_output << rp_obj.index1 << "\t" << rp_obj.index2 << "\t" << rp_obj.start << "\t" << rp_obj.end << "\t";
        }
        else{
            false_output << rp_obj.index1 << "\t" << rp_obj.index2 << "\t" << rp_obj.start << "\t" << rp_obj.end << "\t";
        }
        if(!std::getline(hap_input, line_r)) break;
        IBDSegment rp_obj_new = IBDSegment(line_r);
        if (!(gt_obj == rp_obj_new)) {
            ibd_r.clear();
        }
        rp_obj = rp_obj_new;
    }
    hap_input.close();
    gt_input.close();
    true_output.close();
    false_output.close();
}

std::vector<std::pair<int, int>> overlappingWindows(std::vector<float> cm, float min_seed, int min_markers, int n_threads){
    std::vector<std::pair<int, int>> window_vec;
    float L = cm.back() - cm.front();
    float window_size = ((L - min_seed) / n_threads + min_seed);
    float dist = (L - min_seed) / n_threads;
    float start_gen_pos = 0.0f;
    int i = 0;
    while(i < n_threads){
        float end_gen_pos = start_gen_pos + window_size;
        int start = findInsertionIndex(cm, start_gen_pos);
        int end = findInsertionIndex(cm, end_gen_pos);
        start_gen_pos += dist;
        window_vec.push_back(std::pair<int, int>(start, end));
        i++;
    }
    return window_vec;
}
int minSites(std::vector<float> &cm_mapping, float min_seed){
	int min = cm_mapping.size();
	for(int i = 0; i < cm_mapping.size(); i++){
		int j = i + 1;
		
		while ((j < cm_mapping.size()) && (cm_mapping[j] - cm_mapping[i] < min_seed)){
			j++;
		}
		if((cm_mapping.back() - cm_mapping[i] > min_seed)){
			int n = j - i + 1;
			if (n < min){
				min = n;
			}
		}
		
	}
	return (min == cm_mapping.size()) ? -1 : min;
}