/**
 * Author: Joseph Valencia
 * Date: 12/14/17
 * Bioinformatics Toolsmith Laboratory, University of Tulsa
 * */
#include <string>
#include <algorithm>
#include <vector>
#include <iostream>
#include <limits.h>
#include <fstream>
#include <string.h>
#include <cmath>
#include "GlobAlign.h"

using namespace std;

GlobAlign::GlobAlign(const char * seq1In, int start1In, int end1In, const char * seq2In,
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

void GlobAlign::findAlignment(){

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
            int ans = max(xgapBegin,xgapCont);

            // 0 represents matches/mismatches, 1 represents lowerGap, 2 represents upperGap
            if(ans == xgapBegin){
                lowerDir[at(i,j)] = 0; //backtracking table
            }
            else{
                lowerDir[at(i,j)] =1;
            }

            lowerGap[at(i,j)]= ans; //score table for gap in sequence 2

            //compute values for upperGap
            int ygapBegin = matches[at(i,j-1)]- (gapOpen+gapContinue);
            int ygapCont = upperGap[at(i,j-1)] -gapContinue;

            ans = max(ygapBegin,ygapCont);

            if(ans == ygapBegin){
                upperDir[at(i,j)] = 0;
            }
            else{
                upperDir[at(i,j)] = 2;
            }

            upperGap[at(i,j)]= ans; // score table for gap in seq1
            //compute values for matches
            int matchScore = (seq1[start1+i-1] == seq2[start2+j-1]) ? match : mismatch;

            int matched = matches[at(i-1,j-1)] + matchScore;
            int xgapEnd = lowerGap[at(i-1,j-1)] + matchScore;
            int ygapEnd = upperGap[at(i-1,j-1)] + matchScore;

            ans = max(max(matched,xgapEnd),ygapEnd);

            if(ans == matched){
                matchDir[at(i,j)] = 0;
            }
            else if( ans == xgapEnd){
                matchDir[at(i,j)] =1;
            }
            else{
                matchDir[at(i,j)] =2;
            }

            matches[at(i,j)] = ans; //score table for match/mismatch in seq1 and seq2

        }
    }

    alignmentScore = max(max(matches[at(len1-1,len2-1)], lowerGap[at(len1-1,len2-1)]), upperGap[at(len1-1,len2-1)]); //Score found
    //cout <<"ALIGNMENT SCORE FOUND:"<<alignmentScore<<endl;
    //Begin backtracking code
    int table;
    if(alignmentScore == matches[at(len1-1,len2-1)]){
            table = 0;
    }
    else if(alignmentScore == lowerGap[at(len1-1,len2-1)]){
            table = 1;
    }
    else{
        table = 2;
    }

    vector<char> top;
    vector<char>bottom;

    int i = len1-1;
    int j = len2-1;

    int curr;

    while(i !=0 && j!=0){

        if(table==0){
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

    }

    while(j!=0){
        top.push_back('_');
        bottom.push_back(seq2[j-1]);
        j-=1;
    }

    while(i!=0){
        top.push_back(seq1[i-1]);
        bottom.push_back('_');
        i-=1;
    }


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

int GlobAlign::getScore(){
    return alignmentScore;
}
int GlobAlign::getLength(){
    return bottomString.size();
}
void GlobAlign::printAlignment(){

    /*for(int i = 0;i<topString.size();i++){
        cout << topString[i];
    }*/
    cout <<endl;
    for(int i = 0;i<bottomString.size();i++){
        cout << bottomString[i];
    }
    cout <<endl;

}


double GlobAlign::getIdentity() const
{
    double matches = 0;
    for (int i = 0; i < topString.size(); i++) {
        if (topString[i] == bottomString[i]) {
            matches++;
        }
    }

    cout <<"MATCHES"<<matches;
    return matches / topString.size();
};

 int main(int argc, char *argv[]){

    ifstream ifs;
    
      ifs.open (argv[1], ifstream::in);
      cout<<"FILE OPENED"<<endl;
      char c = ifs.get();
 
      if(c == '>'){
         
         while(c!='\n'){
             c = ifs.get();
             
         }
      }
     
     string string1  ="";
     
      while (ifs.good()) {
        
        
        if (c!='\n'){
        string1+=c;
        }
        c = ifs.get();
      }
    
      ifs.close();
 
 
     ifstream ifs2;
     
     ifs2.open (argv[2], ifstream::in);
     
     c = ifs2.get();
  
     if(c == '>'){
  
         while(c!='\n'){
              c = ifs2.get();
         }
     }
  
     string string2  ="";
      
     while (ifs2.good()) {
         
         if(c!='\n'){
         string2+=c;
         }
         c = ifs2.get();
     }
     
     ifs2.close();
    
      
     
  
     std::transform(string1.begin(),string1.end(),string1.begin(),::toupper);
     std::transform(string2.begin(),string2.end(),string2.begin(),::toupper);
 
 
    //cout <<string1<<endl;
    //cout <<string2<<endl;
    
 
     GlobAlign * align = new GlobAlign(string1.c_str(), 0, string1.size()-1,string2.c_str(), 0,string2.size()-1,1,-1,4,1);
     cout << "FINAL SCORE ="<<align->getScore()<<endl;
     cout << "ALIGNMENT LENGTH =" <<align->getLength()<<endl;
     cout << "IDENTITY"<<align->getIdentity()<<endl;
     }
