/*
 * LCS.cpp
 *
 *  Created on: Dec 4, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "LCS.h"

#include <string>
#include "../exception/InvalidStateException.h"

#include <vector>
#include <iostream>
//test
using namespace std;
using namespace utility;
using namespace exception;

const string LCS::SAME = "same";
const string LCS::OPPOSITE = "opposite";

LCS::LCS(const char * seq1In, int start1In, int end1In, const char * seq2In,
		int start2In, int end2In) {
	seq1 = seq1In;
	start1 = start1In;
	end1 = end1In;

	seq2 = seq2In;
	start2 = start2In;
	end2 = end2In;

	len1 = end1 - start1 + 2;
	len2 = end2 - start2 + 2;
	lenTotal = len1 * len2;

	cTable = new int[lenTotal];
	for (int i = 0; i < lenTotal; i++) {
		cTable[i] = 0;
	}

	bTable = new int[lenTotal];
	for (int i = 0; i < lenTotal; i++) {
		bTable[i] = 0;
	}

	findLcs();
}

LCS::~LCS() {
	delete[] cTable;
	delete[] bTable;
}

void LCS::findLcs() {
	for (int i = 1; i < len1; i++) {
		char base1 = seq1[start1 + i - 1];
		int iM1Index = (i - 1) * len2;
		int iIndex = i * len2;

		for (int j = 1; j < len2; j++) {
			int ijIndex = iIndex + j;

			if (base1 == seq2[start2 + j - 1]) {
				cTable[ijIndex] = cTable[iM1Index + j - 1] + 1;
				bTable[ijIndex] = DIAGONAL;
			} else {
				if (cTable[iM1Index + j] > cTable[iIndex + j - 1]) {
					cTable[ijIndex] = cTable[iM1Index + j];
					bTable[ijIndex] = UP;
				} else {
					cTable[ijIndex] = cTable[iIndex + j - 1];
					bTable[ijIndex] = LEFT;
				}
			}
		}
	}

	// Testing
	/*
	for (int i = 0; i < len1; i++) {
		int iIndex = i * len2;
		for (int j = 0; j < len2; j++) {
			cout << cTable[iIndex + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;

	for (int i = 0; i < len1; i++) {
		int iIndex = i * len2;
		for (int j = 0; j < len2; j++) {
			cout << bTable[iIndex + j] << " ";
		}
		cout << endl;
	}
	cout << endl;
	*/
	// End
}

int LCS::getLenCS(){
	return cTable[lenTotal-1];
}

void LCS::printLcs() {
	vector<char> * rev = new vector<char>();

	int i = len1 - 1;
	int j = len2 - 1;

	while (i != 0 && j != 0) {
		int iIndex = i * len2;
		switch (bTable[iIndex + j]) {
		case DIAGONAL:
			// Test start
			cout << start1 + i - 1 << endl;
			// Test end
			rev->push_back(seq1[start1 + i - 1]);
			i--;
			j--;
			break;
		case UP:
			i--;
			break;
		case LEFT:
			j--;
			break;
		default:
			string msg = "Invalid direction in the bTable: ";
			msg.append(1, bTable[iIndex + j]);
			msg.append(1, '.');
			throw InvalidStateException(msg);
			break;
		}
	}

	int size = rev->size();
	cout << "Rev size is: " << size << endl;

	for (int i = size - 1; i >= 0; i--) {
		cout << /*(int)*/ rev->at(i);
	}
	cout << endl;

	rev->clear();
	delete rev;
}
