#include <segment.hpp>


bool IBDSegment::operator<(const IBDSegment& other){
    if(this->index1 < other.index1){
        return true;
    }
    if(this->index1 == other.index1){
        if (this->index2 <= other.index2){
            return this->index2 >= other.index2;
        }
    }
    return false;
}

bool IBDSegment::operator>(const IBDSegment& other){
    if(this->index1 > other.index1){
        return true;
    }
    if(this->index1 == other.index1){
        if(this->index2 > other.index2){
            return true;
        }
    }
    return false;
}

bool IBDSegment::operator<=(const IBDSegment& other){
    if (this->index1 < other.index1){
        return true;
    }
    if (this->index1 == other.index1){
        return this->index2 <= other.index2;
    }
    return false;
}

bool IBDSegment::operator==(const IBDSegment& other){
    return (this->index1 == other.index1 && this->index2 == other.index2);
}

double IBDSegment::getCoverage(IBDSegment& other){
    // returns propotion of this covered by other
    double prop = static_cast<double>((std::max(0, std::min(this->segment_interval.end, other.segment_interval.end) - std::max(this->segment_interval.start, other.segment_interval.start)))/ ((other.segment_interval.end - other.segment_interval.start)));
    return prop;
}

double getCoverage(int start1, int end1, int start2, int end2){
    double prop = static_cast<double>((std::max(0, std::min(end1, end2) - std::max(start1, start2)))/ ((end2 - start2)));
    return prop;

};