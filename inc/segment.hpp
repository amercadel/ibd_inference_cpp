#pragma once
#include "utils.hpp"

class IBDSegment{
    public:
        int index1;
        int index2;
        int start;
        int end;
        interval segment_interval;


        IBDSegment(int index1, int index2, int start, int end){
            this->index1 = index1;
            this->index2 = index2;
            this->start = start;
            this->end = end;
            this->segment_interval.start = start;
            this->segment_interval.end = end;
        }
        IBDSegment(const std::string& str) {
            int index1 = -1;
            int index2 = -1;
            int start = -1;
            int end = -1;
            if(str.empty() || str == ""){
                this->index1 = index1;
                this->index2 = index2;
                this->start = start;
                this->end = end;
                this->segment_interval.start = start;
                this->segment_interval.end = end;
            }
            else{
                std::vector<std::string> vals = split(str, '\t');
                this->index1 = std::stoi(vals[0]);
                this->index2 = std::stoi(vals[1]);
                this->start = std::stod(vals[2]);
                this->end = std::stod(vals[3]);
                this->segment_interval.start = this->start;
                this->segment_interval.end = this->end;
            }
        }

        bool operator<(const IBDSegment& other);

        bool operator>(const IBDSegment& other);

        bool operator<=(const IBDSegment& other);

        bool operator==(const IBDSegment& other);

        double getCoverage(IBDSegment& other);
};

double getCoverage(int start1, int end1, int start2, int end2);