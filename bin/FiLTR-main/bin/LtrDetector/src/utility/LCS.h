/*
 * LCS.h
 *
 *  Created on: Dec 4, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef LCS_H_
#define LCS_H_

#include <string>

using namespace std;

namespace utility{
class LCS {
private:
	const char * seq1;
	int start1;
	int end1;
	const char * seq2;
	int start2;
	int end2;

	int len1;
	int len2;
	int lenTotal;
	// int lenCS;

	int * cTable;
	int * bTable;
	static const int UP = 1;
	static const int DIAGONAL = 2;
	static const int LEFT = 3;

public:
	static const string SAME;
	static const string OPPOSITE;

	LCS(const char *, int, int, const char *, int, int);
	virtual ~LCS();
	void findLcs();
	int getLenCS();
	void printLcs();
};
}

#endif /* LCS_H_ */
