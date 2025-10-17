/**
 * Author: Joseph Valencia 
 * Date: 12/14/17
 * Bioinformatics Toolsmith Laboratory, University of Tulsa
 * */
#include <string>
#include "../exception/InvalidStateException.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <limits.h>
#include <string.h>
#include <cmath>
#include "GlobAlignE.h"

using namespace std;
using namespace utility;
using namespace exception;

GlobAlignE::GlobAlignE(const char * seq1In, int start1In, int end1In, const char * seq2In,
        int start2In, int end2In, int matchIn, int mismatchIn, int gapOpenIn, int gapContinueIn){

    seq1 = seq1In;
	start1 = start1In;
	end1 = end1In;

    seq2 = seq2In;
	start2 = start2In;
	end2 = end2In;

	len1 = end1 - start1 + 2;
    len2 = end2 - start2 + 2;
    match = matchIn;
    mismatch = mismatchIn;
    gapOpen = gapOpenIn;
    gapContinue = gapContinueIn;
    findAlignment();

}

int GlobAlignE::findAlignment(){

    int shorter = min(len2,len1)-1;
    int lenDiff = abs(len2-len1);
    int maxDiff=0;
    
    if (lenDiff >=1){
        maxDiff += -gapOpen- (lenDiff*gapContinue);
    }
    
    maxDiff+= (mismatch* shorter)-1;

    const int negativeInf = maxDiff; 

    int matches[len1];
    int upperGap[len1];
    int lowerGap[len2];
   
    matches[0]= 0;
    upperGap[0] = negativeInf;
    lowerGap[0] = negativeInf;

    for (int i = 1; i<len1;i++){ 
        upperGap[i] = negativeInf;
        matches[i] = negativeInf;
    }

    for (int j = 1;j<len2;j++){
        lowerGap[j] = negativeInf;
    }

    for( int i = 1;i<len1;i++){

        cout << "matches:";
        for ( int x =0;x<len1;x++){
            cout <<matches[x];
        }
        cout <<endl;
        cout << "upper:";
        for ( int x =0;x<len1;x++){
            cout <<upperGap[x];
        }
        cout << endl;
        cout << "lower";
        for(int y = 0;y<len2;y++){
            cout <<lowerGap[y];
        }
        cout <<endl;

        
        for(int j = 1;j<len2;j++){
       
        //compute values for upperGap
        int ygapBegin = matches[i] -(gapOpen+gapContinue);
        int ygapCont = upperGap[i] -gapContinue;
        cout << "yGapBegin: " <<ygapBegin<<endl;
        cout << "yGapCont: " << ygapCont<<endl;

        int ans = max(ygapBegin,ygapCont);
        upperGap[i] = ans;


         //compute values for match/mismatch
         int matchScore = (seq1[start1+i-1] == seq2[start2+j-1]) ? match : mismatch;
         
         int matched = matches[i-1] + matchScore;
         int xgapEnd = lowerGap[j-1] + matchScore;
         int ygapEnd = upperGap[i-1] + matchScore;
         cout << "matched " <<matched<<endl;
         cout << "xgapEnd: " << xgapEnd<<endl;
         cout << "ygapEnd: " << ygapEnd<<endl;
         ans = max(max(matched,xgapEnd),ygapEnd);
         matches[i] =ans;
 

        //compute values for lowerGap
        int xgapBegin = matches[j-1] -(gapOpen+gapContinue);
        int xgapCont = lowerGap[j]- gapContinue;
        cout << "xgapBegin: " << xgapBegin<<endl;
        cout << "xgapCont: " << xgapCont<<endl;
        ans = max(xgapBegin,xgapCont); 
        lowerGap[j] = ans;

        }

    }
    return max(max(matches[len1-1], lowerGap[len2-1]), upperGap[len1-1]);
                   

}

int main(){
    const char * string1 = "T";
    const char * string2 = "TAG";

    string str1 = string(string1);
    string str2 = string(string2);
    
    GlobAlignE * align = new GlobAlignE(string1, 0, str1.size(),string2, 0,str2.size(),1,-1,4,1);
    cout << align->findAlignment()<<endl;
}
