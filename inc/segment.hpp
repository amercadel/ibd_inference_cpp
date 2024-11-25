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
                if(vals[0].at(0) == 't'){
                    this->index1 = (std::stoi(split(vals[0], '_')[1]) * 2) + (std::stoi(vals[1]) - 1);
                    this->index2 = (std::stoi(split(vals[2], '_')[1]) * 2) + (std::stoi(vals[3]) - 1);
                    this->start = std::stoi(vals[5]);
                    this->end = std::stoi(vals[6]);
                    this->segment_interval.start = std::stoi(vals[5]);
                    this->segment_interval.end = std::stoi(vals[6]);
                }
                else{
                    this->index1 = std::stoi(vals[0]);
                    this->index2 = std::stoi(vals[1]);
                    this->start = std::stoi(vals[2]);
                    this->end = std::stoi(vals[3]);
                    this->segment_interval.start = this->start;
                    this->segment_interval.end = this->end;
                }
                
            }
        }

        bool operator<(const IBDSegment& other);

        bool operator>(const IBDSegment& other);

        bool operator<=(const IBDSegment& other);

        bool operator==(const IBDSegment& other);

        bool operator!=(const IBDSegment& other);

        float getCoverage(const IBDSegment& other);

        void display();
};

float getCoverage(int start1, int end1, int start2, int end2);