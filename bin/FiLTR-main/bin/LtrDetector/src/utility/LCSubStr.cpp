/*
 * LCSubStr.cpp
 *
 *  Created on: Dec 19, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "LCSubStr.h"
#include "Util.h"
#include "Location.h"

// Delete start
#include <iostream>
// Delete end

namespace utility {

LCSubStr::LCSubStr(const char * seq1In, ILocation * loc1, const char * seq2In,
		ILocation * loc2) {

	int start1In = loc1->getStart();
	int end1In = loc1->getEnd();
	int start2In = loc2->getStart();
	int end2In = loc2->getEnd();

	if (end1In - start1In > end2In - start1In) {
		seq1 = seq1In;
		start1 = start1In;
		end1 = end1In;
		len1 = end1 - start1 + 2;

		seq2 = seq2In;
		start2 = start2In;
		end2 = end2In;
		len2 = end2 - start2 + 2;

		isFirstShorter = false;
	} else {
		seq1 = seq2In;
		start1 = start2In;
		end1 = end2In;
		len1 = end1 - start1 + 2;

		seq2 = seq1In;
		start2 = start1In;
		end2 = end1In;
		len2 = end2 - start2 + 2;

		isFirstShorter = true;
	}

	lenTotal = 2 * len2;
	cTable = new int[lenTotal];

	for (int i = 0; i < lenTotal; i++) {
		cTable[i] = 0;
	}
	lenCS = 0;

	subStr1 = new vector<ILocation *>();
	subStr2 = new vector<ILocation *>();

	findLCSubStr();
}

LCSubStr::~LCSubStr() {
	delete[] cTable;

	Util::deleteInVector(subStr1);
	delete subStr1;

	Util::deleteInVector(subStr2);
	delete subStr2;
}

void LCSubStr::findLCSubStr() {
	int iM1Index = 0;
	int iIndex = len2;

	for (int i = 1; i < len1; i++) {
		char base1 = seq1[start1 + i - 1];

		for (int j = 1; j < len2; j++) {
			int ijIndex = iIndex + j;
			if (base1 == seq2[start2 + j - 1]) {
				cTable[ijIndex] = cTable[iM1Index + j - 1] + 1;

				if (cTable[ijIndex] == lenCS) {
					subStr1->push_back(
							new Location(start1 + i - lenCS, start1 + i - 1));
					subStr2->push_back(
							new Location(start2 + j - lenCS, start2 + j - 1));
				} else if (cTable[ijIndex] > lenCS) {
					lenCS = cTable[ijIndex];

					Util::deleteInVector(subStr1);
					Util::deleteInVector(subStr2);

					subStr1->push_back(
							new Location(start1 + i - lenCS, start1 + i - 1));
					subStr2->push_back(
							new Location(start2 + j - lenCS, start2 + j - 1));
				}
			} else {
				cTable[ijIndex] = 0;
			}
		}

		if (i != len1 - 1) {
			for (int h = 0; h < len2; h++) {
				cTable[h] = cTable[len2 + h];
			}
		}
	}

// Test start
	// for (int i = 0; i < subStr1->size(); i++) {
	//	cout << "STD 1: " << subStr1->at(i)->toString() << endl;
	//	cout << "STD 2: " << subStr2->at(i)->toString() << endl;
	// }

	// for (int i = 0; i < lenTotal; i++) {
	//	cout << cTable[i];
	// }
	// cout << endl;
// Test end

}

int LCSubStr::getLenCS() {
	return lenCS;
}

vector<vector<ILocation *> *> * LCSubStr::getCSubStr() {
	vector<vector<ILocation *> *> * r = new vector<vector<ILocation *> *>();
	if (isFirstShorter) {
		r->push_back(subStr2);
		r->push_back(subStr1);
	} else {
		r->push_back(subStr1);
		r->push_back(subStr2);
	}

	return r;
}

}
