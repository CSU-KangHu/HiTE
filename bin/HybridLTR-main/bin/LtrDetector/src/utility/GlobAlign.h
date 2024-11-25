
/**
 * Author: Joseph Valencia
 * Date: 12/14/17
 * Bioinformatics Toolsmith Laboratory, University of Tulsa
 * */
#ifndef Glob_Align_H_
#include <string>

using namespace std;

class GlobAlign{

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
    int alignmentScore;
    string topString;
    string bottomString;
public:
    GlobAlign(const char*,int,int,const char *,int,int, int,int,int,int);
   // virtual  GlobAlign();
    double getIdentity() const;
    void findAlignment();
    void printAlignment(); //display LocAlign
    int getScore();
    int getLength();

};
#endif
