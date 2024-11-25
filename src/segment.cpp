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

float IBDSegment::getCoverage(const IBDSegment& other){
    // returns proportion of this covered by other
    int overlap_start = std::max(this->start, other.start);
    int overlap_end = std::min(this->end, other.end);
    if (overlap_start >= overlap_end) {
        return 0.0f;
    }
    float overlap_length = static_cast<float>(overlap_end - overlap_start);
    float this_length = static_cast<float>(this->end - this->start);
    return overlap_length / this_length;
}

float getCoverage(int start1, int end1, int start2, int end2){
    float prop = static_cast<float>((std::max(0, std::min(end1, end2) - std::max(start1, start2)))/ ((end2 - start2)));
    return prop;

};

void IBDSegment::display(){
    std::cout << this->index1 << "\t" << this->index2 << "\t" << this->start << "\t" << this->end << std::endl;
}