#ifndef Loc_Align_H_
#include <string>
#include <tuple>

using namespace std;

namespace utility{

class LocAlign{

private:
    const char * seq1; //first sequence to be aligned
    int start1;
    int end1;
    const char * seq2;//second sequence to be aligned
    int start2;
    int end2;
    int len1;
    int len2;
    int lenTotal;
    int match; //score for base pair match
    int mismatch;//score for base pair mismatch
    int gapOpen; //cost to open a gap
    int gapContinue; //cost to continue a gap
    int * matches;
    int * upperGap;
    int * lowerGap;
    int * matchLen;
    int * upperLen;
    int * lowerLen;
    int * matchId;
    int * upperId;
    int * lowerId;
    int alignmentScore;
    int alignmentLength;
    int totalMatches;
    int queryStart;
    int queryEnd;
    int referenceStart;
    int referenceEnd;
    string topString; 
    string bottomString;
public:
    LocAlign(const char*,int,int,const char *,int,int, int,int,int,int);
    
    //virtual ~LocAlign();
    void findAlignment();
    tuple<int,int,int>findMax(int *);
    double getIdentity();
    int getLength();
    void printAlignment(); //display LocAlign
    int getScore(); 
    int getLengthAlignment();
    int getQueryStart();
    int getQueryEnd();
    int getReferenceStart();
    int getReferenceEnd();
  
};
}
#endif 