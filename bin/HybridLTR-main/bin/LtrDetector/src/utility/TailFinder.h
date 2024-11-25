/*
 * TailFinder.h
 *
 *  Created on: Nov 27, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef TAILFINDER_H_
#define TAILFINDER_H_

#include <vector>
#include <string>

#include "ILocation.h"

using namespace std;

namespace utility {
class TailFinder {
private:
	const string * seq;
	ILocation * loc;
	int whichTail;
	int win;
	int minLen;

	int seedLen;
	int gapLen;

	vector<int> * tail;

	void findMark();
	void findMarkA(string *, vector<int> *, int, int);
	void findMarkP(string *, vector<int> *, int, int);


public:
	string prettyFormatChrom(string *);
	static const int MARK_A = 1;
	static const int MARK_P = 2;
	TailFinder(const string *, ILocation *, int, int,int,int, int);
	virtual ~TailFinder();
	vector<int> * getTail();
	bool isTailFound();
};
}

#endif /* TAILFINDER_H_ */
