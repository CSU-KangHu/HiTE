/*
 * ScorerTr.h
 *
 *  Created on: Nov 30, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef SCORERTR_H_
#define SCORERTR_H_

#include "ForwardTr.h"
#include "../nonltr/KmerHashTable.h"
#include "../nonltr/ChromosomeOneDigit.h"

#include <vector>

using namespace nonltr;

namespace tr {
class ScorerTr {
private:
	KmerHashTable <int, int> * kmerTable; 
	
	static const int INITIAL_VALUE;
	static const int INITIAL_SCORE;

	ChromosomeOneDigit *chrom;
	int k;
	int min;
	int max;
//	int minPlateauLen;
   // int diffThresh;
  //  int gapTol;
	void medianSmooth();
	int findMedian(int,int);
	void score();
	void scoreNew();
	
	vector<int> * scoreList;
	std::string csvFileName;
	

public:
	ScorerTr(ChromosomeOneDigit *, int, int, int);
	virtual ~ScorerTr();
	vector<int>* getScores();
	int getInitialScore();
	//void outputScores(int, int);
	void scoresFormat(int,int);
	//void bedFormat(int,int);
	//vector<std::tuple<char,int,int>> getSpikes();
};
}

#endif /* SCORERTR_H_ */
