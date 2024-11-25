/*
  Created by Joseph Valencia 21 February 2018
*/
#ifndef MATCHTR_H_
#define MATCHTR_H_

#include <vector>
#include <tuple>
#include "BackwardTr.h"
#include "Candidate.h"
//using namespace nonltr;
using namespace std;
namespace tr{

class MatchTr{

private:
    
    int initValue;
    int k;
    string bedFileName;
    int min;
    int max;
    int ltrMin;
    int minPlateauLen;
    int diffThresh;
    int gapTol;
    int identity;

    bool isMatch(int,int);
    void forwardMerge();
	void backwardMerge();
	void cleanAndMerge();
    void medianSmooth();

    vector<int> * scoreList;
	vector<std::tuple<char,int,int>> spikes;
	vector<Candidate *>  plateaus;
	vector<BackwardTr *> * bList;

    void smoothAndCollect();
	int findMedian(int,int);


public:
    vector <BackwardTr*> * getRepeatCandidates();
    MatchTr(vector<int> *,int,int,int,int,int,int,int,int);
    virtual ~MatchTr();
    void bedFormat(int,int);
    void printFinalScores(int, int);
    void adjustBounds();
    vector<int> * getScoreList();
};

}
#endif 
