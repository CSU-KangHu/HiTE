#include <iostream>
#include <cstring>
#include <string>

#include <cstdlib>
#include "../utility/GlobAlignE.h"

using namespace std;
using namespace utility;

int main(int argc, char * argv []){

    string first = string(argv[1]);
    const char * seq1 = first.c_str();
    int start1 = 0;
    int end1 = first.length();
    string second = string(argv[2]);
    
    const char * seq2 = second.c_str();
    int start2 = 0;
    int end2 = second.length();
    int match = atoi(argv[3]);
    int mismatch = atoi(argv[4]);
    int gapOpen = atoi(argv[5]);
    cout<<gapOpen<<endl;
    int gapContinue = atoi(argv[6]);
   

    GlobAlignE * aligner = new GlobAlignE(seq1,start1,end1,seq2,start2,end2,match,mismatch,gapOpen,gapContinue);
    
    cout<<aligner->getIdentity()<<endl;
   // delete aligner;
    return 0;

}