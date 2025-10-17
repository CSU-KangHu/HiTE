/*
 * FilterTr.h
 *
 *  Created on: Dec 14, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef FILTERTR_H_
#define FILTERTR_H_

#include "BackwardTr.h"
#include "LtrTe.h"

#include <string>
#include <vector>

using namespace std;

namespace tr {

class FilterTr {
private:
	const string * seq;
	const char * cSeq;
	vector<BackwardTr *> * bList;
	vector<LtrTe *> * teList;
	string bedFileName;
	string name;
	int k;
	int init;
	int ltrLen;
	int ltrId;
	int ltrSep;
	int tsdW;
	int tsdT;
	// int tailW;
	int min;
	int max;
	int maxLtrLen;
	int minLtrLen;

	const int tailT = 12;

	bool canUseLtr;
	bool canUseSine;
	bool canUsePpt;
	bool canUseTsd;
	bool canUseLength;
	bool canUseDNA;
	
    void tightenBounds();
	void adjust();
	void filter();
	void fillTeList();
	void filterAcc2Ltr();
	void filterAcc2Sine();
	void filterAcc2Tsd();
	void filterAcc2Ppt();
	void filterAcc2Length();
	void filterAcc2DNA();

	int calculateTailWindow(double,int,int);
	void removeOverlaps();

public:
	FilterTr(string,const string *, vector<BackwardTr *> *, int,int,int,int,int,int);
	virtual ~FilterTr();
	vector<LtrTe *> * getTeList();
	bool orderFunction(BackwardTr *, BackwardTr *);
	void bedFormat(int,int);
	void fullFormat(int,int);
	string convertNucleotides(string);

	
};

} /* namespace tr */
#endif /* FILTERTR_H_ */
