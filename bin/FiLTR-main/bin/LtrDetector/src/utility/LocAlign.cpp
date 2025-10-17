/**
 * Author: Joseph Valencia
 * Date: 12/14/17
 * Bioinformatics Toolsmith Laboratory, University of Tulsa
 * */
#include <string>
#include <algorithm>
#include <vector>
#include <tuple>
#include <iostream>
#include <limits.h>
#include <fstream>
#include <string.h>
#include <cmath>
#include "LocAlign.h"

using namespace std;
using namespace utility;



LocAlign::LocAlign(const char *seq1In, int start1In, int end1In, const char *seq2In,
                   int start2In, int end2In, int matchIn, int mismatchIn, int gapOpenIn, int gapContinueIn)
{

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

void LocAlign::findAlignment(){

    int shorter = min(len2,len1) -1;
    int lenDiff = abs(len2-len1);
    int maxDiff=0;

    if (lenDiff >=1){
        maxDiff += -gapOpen- (lenDiff*gapContinue);
    }

    maxDiff+= (mismatch* shorter)-1;

    const int negativeInf = maxDiff; // negative infinity = one less than minimum score given the length of two strings

    auto at = [&](int a, int b) -> int { return a * len2 + b; };

    int *matches = new int[len1*len2];
    int *lowerGap= new int[len1*len2];
    int *upperGap= new int[len1*len2];
    int *matchDir= new int[len1*len2];
    int *lowerDir= new int[len1*len2];
    int *upperDir= new int[len1*len2];

   

    matches[at(0,0)] = 0;
    lowerGap[at(0,0)] = negativeInf;
    upperGap[at(0,0)] = negativeInf;
    matchDir[at(0,0)] =0;
    lowerDir[at(0,0)] =0;
    upperDir[at(0,0)] =0;

    

    //initialize first col
    for (int i = 1; i<len1;i++){

        lowerGap[at(i,0)] = (-gapOpen)- (i*gapContinue);
        upperGap[at(i,0)] = negativeInf;
        matches[at(i,0)] = negativeInf;

        matchDir[at(i,0)] =0;
        lowerDir[at(i,0)] =0;
        upperDir[at(i,0)] =0;
    }

    //initialize first row
    for( int j = 1;j<len2;j++){

        upperGap[at(0,j)] = (-gapOpen)-(j*gapContinue);
        lowerGap[at(0,j)] = negativeInf;
        matches[at(0,j)] = negativeInf;


        matchDir[at(0,j)] =0;
        lowerDir[at(0,j)] =0;
        upperDir[at(0,j)] =0;
    }
 
 
    for( int i = 1;i<len1;i++){


        for(int j = 1;j<len2;j++){


            //compute values for lowerGap
            int xgapBegin = matches[at(i-1,j)] -(gapOpen+gapContinue);
            int xgapCont = lowerGap[at(i-1,j)]- gapContinue;
            int ans = max(0,max(xgapBegin,xgapCont));

            // 0 represents from zero, 1 represents lowerGap, 2 represents upperGap ,3 represents match/mismatches
            if(ans == xgapBegin){
                lowerDir[at(i,j)] = 3; //backtracking table
            }
            else if(ans ==xgapCont){
                lowerDir[at(i,j)] =1;
            }
            else{
                lowerDir[at(i,j)] = 0;
            }

            lowerGap[at(i,j)]= ans; //score table for gap in sequence 2

            //compute values for upperGap
            int ygapBegin = matches[at(i,j-1)]- (gapOpen+gapContinue);
            int ygapCont = upperGap[at(i,j-1)] -gapContinue;

            ans = max(0,max(ygapBegin,ygapCont));

            if(ans == ygapBegin){
                upperDir[at(i,j)] = 3;
            }

            else if (ans == ygapCont)
                {
                    upperDir[at(i, j)] = 2;
                }
            else{
                upperDir[at(i,j)] = 0;
            }

            upperGap[at(i,j)]= ans; // score table for gap in seq1
            //compute values for matches
            int matchScore = (seq1[start1+i-1] == seq2[start2+j-1]) ? match : mismatch;

            int matched = matches[at(i-1,j-1)] + matchScore;
            int xgapEnd = lowerGap[at(i-1,j-1)] + matchScore;
            int ygapEnd = upperGap[at(i-1,j-1)] + matchScore;

            ans = max(0, max(max(matched,xgapEnd),ygapEnd));

            if(ans == matched){
                matchDir[at(i,j)] = 3;
            }
            else if( ans == xgapEnd){
                matchDir[at(i,j)] =1;
            }
            else if (ans ==ygapEnd){
                matchDir[at(i,j)] =2;
            }
            else{
                matchDir[at(i,j)] = 0;
            }

            matches[at(i,j)] = ans; //score table for match/mismatch in seq1 and seq2

        }
    }



    int maxMatchi;
    int maxMatchj;
    int maxMatchVal;

    int maxUpperi;
    int maxUpperj;
    int maxUpperVal;

    int maxLoweri;
    int maxLowerj;
    int maxLowerVal;

    std::tie(maxMatchVal,maxMatchi,maxMatchj) = findMax(matches);
    std::tie(maxUpperVal,maxUpperi,maxUpperj) = findMax(upperGap);
    std::tie(maxLowerVal,maxLoweri,maxLowerj) = findMax(lowerGap);


    alignmentScore = max(max(maxMatchVal,maxUpperVal),maxLowerVal); //Score found
    //cout <<"ALIGNMENT SCORE FOUND:"<<alignmentScore<<endl;
    //Begin backtracking code

    int table;

    int soli;
    int solj;

    if(alignmentScore == maxMatchVal){
            table = 3;
            soli = maxMatchi;
            solj = maxMatchj;
    }
    else if(alignmentScore == maxLowerVal){
            table = 1;
            soli = maxLoweri;
            solj = maxLowerj;
    }
    else if (alignmentScore == maxUpperVal){
        table = 2;
        soli = maxUpperi;
        solj = maxUpperj;
    }

  

    vector<char> top;
    vector<char>bottom;

    int i = soli;
    int j = solj;

    int curr;

    while(i !=0 && j!=0){

        if(table==3){
            curr = matchDir[at(i,j)];
            top.push_back(seq1[i-1]);
            bottom.push_back(seq2[j-1]);
            i-=1;
            j-=1;
            table = curr;

        }
        else if(table==1){
           curr = lowerDir[at(i,j)];
           top.push_back(seq1[i-1]);
           bottom.push_back('_');
           i-=1;
           table = curr;
        }
        else if(table ==2){
           curr = upperDir[at(i,j)];
           top.push_back('_');
           bottom.push_back(seq2[j-1]);
           j-=1;
           table = curr;
        }
        else if (table ==0){
            i++;
            j++;
            break;
           
        }

    }


    queryStart = i;
    queryEnd = soli;

    referenceStart = j;
    referenceEnd = solj;
  
   // cout<<"Top alignment = "<<i<<":"<<soli<<endl;
    //cout<<"Bottom alignment = "<<j<<":"<<solj<<endl;

    reverse(top.begin(),top.end());
    reverse(bottom.begin(),bottom.end());

    string t ="";
    string b = "";
    for (std::vector<char>::iterator it1=top.begin(); it1!=top.end(); ++it1){
            t+= *it1;
    }


    for (std::vector<char>::iterator it2=bottom.begin(); it2!=bottom.end(); ++it2){
           b+= *it2;
    }

  

    topString = t; //upper part of LocAlign
    bottomString = b;//lower part of LocAlign
    delete[] matches;
    delete[] lowerGap;
    delete[] upperGap;
    delete[] matchDir;
    delete[] lowerDir;
    delete[] upperDir;
}
tuple<int,int,int> LocAlign::findMax( int * array){

    auto at = [&](int a, int b) -> int { return a * len2 + b; };

    int max = array[at(0,0)];
    int maxi =0;
    int maxj =0;

    for (int i =0;i<len1;i++){

        for(int j =0;j<len2;j++){

            int curr = array[at(i,j)];

            if (curr>max){
                max = curr;
                maxi = i;
                maxj = j;
            }
        }
       
    }
    return make_tuple(max,maxi,maxj);

}

int LocAlign::getScore(){
    return alignmentScore;
}
int LocAlign::getLength(){
    return bottomString.size();
}
int LocAlign::getQueryEnd(){
    return queryEnd;
}
int LocAlign::getQueryStart(){
    return queryStart;
}
int LocAlign::getReferenceEnd(){
    return referenceEnd;
}
int LocAlign::getReferenceStart(){
    return referenceStart;
}
void LocAlign::printAlignment(){

    /*for(int i = 0;i<topString.size();i++){
        cout << topString[i];
    }*/
    cout <<endl;
    for(int i = 0;i<bottomString.size();i++){
        cout << bottomString[i];
    }
    cout <<endl;

    cout << endl;
    for (int i = 0; i < topString.size(); i++)
    {
        cout << topString[i];
    }
    cout << endl;
}


double LocAlign::getIdentity() 
{
    double matches = 0;

    for (int i = 0; i < topString.size(); i++) {
        if (topString[i]!='_' && bottomString[i]!='_' && topString[i] == bottomString[i]) {
            matches++;
        }
    }

    //cout <<"MATCHES"<<matches<<"/"<<topString.size()<<endl;

    return matches / topString.size();
};

