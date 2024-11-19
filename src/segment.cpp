#include <segment.hpp>
#include <iostream>


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
bool IBDSegment::operator!=(const IBDSegment& other){
    return (this->index1 != other.index1 || this->index2 != other.index2);
}

float IBDSegment::getCoverage(IBDSegment& other){
    // returns propotion of this covered by other
    float prop = static_cast<float>((std::max(0, std::min(this->segment_interval.end, other.segment_interval.end) - std::max(this->segment_interval.start, other.segment_interval.start)))/ ((other.segment_interval.end - other.segment_interval.start)));
    return prop;
}

float getCoverage(int start1, int end1, int start2, int end2){
    float prop = static_cast<float>((std::max(0, std::min(end1, end2) - std::max(start1, start2)))/ ((end2 - start2)));
    return prop;

};

void IBDSegment::display(){
    std::cout << this->index1 << "\t" << this->index2 << "\t" << this->start << "\t" << this->end << std::endl;
}