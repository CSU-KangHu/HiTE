#include "Candidate.h"
#include <cstdlib>
#include <iostream>
using namespace std;

namespace tr{

    Candidate::Candidate(int startIn,int endIn,int heightIn){

        start = startIn;
        end = endIn;
        height = heightIn;
    }

    Candidate::~Candidate(){

    }

     int Candidate::getStart(){

        return start;
    }

    int Candidate::getEnd(){
        return end;
    }

    int Candidate::getHeight(){
        return height;
    }
    int Candidate::getAbsHeight(){
        return abs(height);
    }

    void Candidate::printCandidate(){
        cout<<"start: "<<getStart()<<" end: "<<getEnd()<<" height: "<<getHeight()<<endl;
    }
}