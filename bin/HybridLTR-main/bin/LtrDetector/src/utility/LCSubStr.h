/*
 * LCSubStr.h
 *
 *  Created on: Dec 19, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef LCSUBSTR_H_
#define LCSUBSTR_H_

#include <vector>

#include "ILocation.h"

using namespace std;

namespace utility {

class LCSubStr {
private:
	const char * seq1;
	int start1;
	int end1;
	int len1;

	const char * seq2;
	int start2;
	int end2;
	int len2;

	bool isFirstShorter;
	int lenTotal;

	int * cTable;
	int lenCS;

	vector< ILocation * > * subStr1;
	vector< ILocation * > * subStr2;

	void findLCSubStr();

public:
	LCSubStr(const char *, ILocation *, const char *, ILocation *);
	virtual ~LCSubStr();
	int getLenCS();
	vector< vector< ILocation * > *> * getCSubStr();
};

}

#endif /* LCSUBSTR_H_ */
