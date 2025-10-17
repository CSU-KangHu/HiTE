/*
 * DetectorTr.h
 *
 *  Created on: Dec 6, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef DETECTORTR_H_
#define DETECTORTR_H_

#include <vector>
#include <string>
#include "BackwardTr.h"
#include "ForwardTr.h"

using namespace std;

namespace tr {

class DetectorTr {
private:
	vector<int>* scoreList;
	int initValue;
	int gap;

	vector<BackwardTr *> * bList;
	vector<ForwardTr *> * fList;

	void findTr();
	void findTrHelper(int, vector<int> *);
	void findMatch(int, int, int, int);
	void sortLastTr();
	void matchBackwardTr(BackwardTr &);
	int checkMatch(int, int);
	void addToBList(BackwardTr *);

public:
	DetectorTr(vector<int>*, int);
	virtual ~DetectorTr();
	vector<BackwardTr *> * getBList();
};

} /* namespace tr */
#endif /* DETECTORTR_H_ */
