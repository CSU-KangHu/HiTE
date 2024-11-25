/*
 * DetectorTr.cpp
 *
 *  Created on: Dec 6, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "DetectorTr.h"
#include "BackwardTr.h"
#include "ForwardTr.h"
#include "../utility/Util.h"
#include "../exception/InvalidStateException.h"
#include "../exception/InvalidInputException.h"

// ToDo: delete the following include statement
#include <iostream>
// end

using namespace std;
using namespace utility;
using namespace exception;

namespace tr {

DetectorTr::DetectorTr(vector<int>* scoreListIn, int initValueIn) {
	scoreList = scoreListIn;
	initValue = initValueIn;
	gap = 100;
	bList = new vector<BackwardTr *>();
	fList = new vector<ForwardTr *>();

	findTr();

	// Post condition:
	// Test if the backward list is correctly sorted and
	// No two elements are overlapping
	for (int j = 1; j < bList->size(); j++) {
		BackwardTr * b1 = bList->at(j);
		BackwardTr * b2 = bList->at(j - 1);

		if (b1->getE1() < b2->getE1()) {
			string msg("List is not increasing! ");
			msg.append("J-1: ");
			msg.append(b2->toString());
			msg.append(" J: ");
			msg.append(b1->toString());
			msg.append(".");
			throw InvalidStateException(msg);
		}

		if (b1->isOverlapping(b2)) {
			string msg("There is overlapping matches in the tr list! ");
			msg.append("J-1: ");
			msg.append(b2->toString());
			msg.append(" J: ");
			msg.append(b1->toString());
			msg.append(".");

			throw InvalidStateException(msg);
		}
	}

	// Test start
	cout << "BList size: " << bList->size() << endl;
	cout << "FList size: " << fList->size() << endl;
	// Test end
}

DetectorTr::~DetectorTr() {
	Util::deleteInVector(bList);
	delete bList;

	Util::deleteInVector(fList);
	delete fList;
}
/*vector<std::tuple<int,int>> DetectorTr::findPlateaus(int diffTolerance){
    
	int len = scoreList->size();
	vector<std::tuple<int,int>> positions;

    int level = initValue;
	int plateauBegin = 0;
	int plateauEnd = 0;


	for(int i =0; i<len;i++){
		
		int currVal = scoreList->at(i);
		
		if(currVal!= initValue){
            
			level = currVal;
			plateauBegin = i;
			plateauEnd = i;
			bool extend = true;
			int gap = 0;

			while(extend){
				
				if (i+1 >=len-1){
					break;
				}
                int nextVal = scoreList->at(i+1);
				int plateauLength = plateauEnd-plateauBegin;

				
				if(abs(nextVal-currVal)>=diffTolerance){
					https://code.visualstudio.com/docs/editor/codebasics
					if(plateauLength < MIN_STRETCH){
                        extend = false;
					}
					else if (gap<=MAX_GAP){
						gap+=1;
					}
					else{
						extend = false;
					}

				}
				else{
                    if(abs(nextVal-currVal) <= diffTolerance){
						i++;
						plateauEnd = i;
						gap=0;
					}
				}
				n != initValue &&n
			}
			if(abs(plateauEnd-plateauBegin)>=PLATEAU_MIN_LEN){
			  std::tuple<int,int> plateauBounds = std::make_tuple(plateauBegin,plateauEnd);
			  positions.push_back(plateauBounds);
			}

			
		} int matchLoc = currStart+scoreList->at(currStart);

				if(matchLoc ==)
		
	}
	return positions;
}
*/

void DetectorTr::findTr() {
	int len = scoreList->size();

	for (int i = 0; i < len; i++) {
		if (scoreList->at(i) != initValue) {
			// Fix: coor is now allocated on the stack not on the heap
			vector<int> coor = vector<int>();
			findTrHelper(i, &coor);
			int start1 = coor.at(0);
			int end1 = coor.at(1);
			int start2 = scoreList->at(i+start1);
			int end2 = scoreList->at(i+end1);
			findMatch(start1, end1, start2, end2);
			i = checkMatch(coor.at(0), coor.at(1));
		}
	}
}

	void DetectorTr::findTrHelper(int i, vector<int> *coor)
	{
		int s = scoreList->at(i);
		if (s == initValue)
		{
			string msg("The score at the index must be initialized.");
			throw InvalidInputException(msg);
		}

		int len = scoreList->size();
		int gapEnd = i + gap - 1;
		if (gapEnd >= len)
		{
			gapEnd = len - 1;
		}
		int start = i;
		int end = i;

		for (int j = i + 1; j <= gapEnd; j++)
		{
			int n = scoreList->at(j);

			if (abs(n - s) <= 5) //DIFF_TOLERANCE
			{
				gapEnd = j + gap - 1;
				if (gapEnd >= len)
				{
					gapEnd = len - 1;
				}
				end = j;
				s = n;
			}
		}
		coor->push_back(start);
		coor->push_back(end);
	}

/*void DetectorTr::findTrHelper(int i, vector<int> * coor) {
	int s = scoreList->at(i);
	if (s == initValue) {
		string msg("The score at the index must be initialized.");
		throw InvalidInputException(msg);
	}

	int len = scoreList->size();
	int gapEnd = i + gap - 1;
	if (gapEnd >= len) {
		gapEnd = len - 1;
	}
	int start = i;
	int end = i;
	for (int j = i + 1; j <= gapEnd; j++) {
		int n = scoreList->at(j);

		if (n != initValue && n >= s && j-1 <= gap) {
			gapEnd = j + gap - 1;
			if (gapEnd >= len) {
				gapEnd = len - 1;
			}
			end = j;
			s = n;
		}
	}
	coor->push_back(start);
	coor->push_back(end);
}*/

void DetectorTr::findMatch(int s1, int e1, int s2, int e2) {
	
	int size = fList->size();
	bool isMatchFound = false;

	if (s1>0) {
		fList->push_back(new ForwardTr(s1, e1, s2, e2));
	} else {
		BackwardTr b(s1, e1, s2, e2);

		for (int i = size - 1; i >= 0; i--) {
			ForwardTr * f = fList->at(i);

			if (b.isOverlapping(f)) {
				BackwardTr hani;
				//b.merge(f, hani);
				matchBackwardTr(hani); /* This call deals with a list of backward TRs*/

				fList->erase(fList->begin() + i);
				isMatchFound = true;
			}
		}

		if (!isMatchFound) {
			matchBackwardTr(b);
		}
	}
	// Update the size of fList
	size = fList->size();

	// Match the forward TRs that will not be matched with new backward TRs.
	// The parameter s1 is the last scanned nucleotide.
	// If the end, the end of the right LTR, of the forward TRs is less than s1,
	// then no backward match will be found for those forward TRs. In this case
	// match those forward TRs with the already-matched backward TRs.
	vector<ForwardTr *> * tempList = new vector<ForwardTr *>();
	for (int g = 0; g < size; g++) {
		ForwardTr * f = fList->at(g);
		if (f->getE2() < s1) {
			BackwardTr hani;
			f->flip(hani);
			matchBackwardTr(hani); /* This call deals with a list of backward TRs*/
		} else {
			tempList->push_back(new ForwardTr(*f));
		}
	}

	// Test start
	/*
	 cout << "Temp size: " << tempList->size() << " Temp capacity: "
	 << tempList->capacity() << endl;
	 cout << "FList size: " << fList->size() << " FList capacity: "
	 << fList->capacity() << endl;
	 */
	// Test end
	// Free memory used by the old list
	Util::deleteInVector(fList);
	delete fList;

	// Assign the new list to the old pointer
	fList = tempList;
}

void DetectorTr::addToBList(BackwardTr * b) {
	bList->push_back(b);
	sortLastTr();
}

void DetectorTr::sortLastTr() {
	int size = bList->size();
	if (size > 0) {
		for (int j = size - 2; j >= 0; j--) {
			BackwardTr * b1 = bList->at(j + 1);
			BackwardTr * b2 = bList->at(j);
			if (b1->getE1() < b2->getE1()) {
				BackwardTr * temp = b2;
				(*bList)[j] = b1;
				(*bList)[j + 1] = temp;
			} else {
				break;
			}
		}
	}
}

/**
 * @@
 * Candidate for memory leak
 */
void DetectorTr::matchBackwardTr(BackwardTr & b) {
	vector<BackwardTr *> overlapList = vector<BackwardTr *>();
	vector<int> eraseList = vector<int>();

	int size = bList->size();
	for (int i = size - 1; i >= 0; i--) {
		BackwardTr * b1 = bList->at(i);
		if (b.isOverlapping(b1)) {
			overlapList.push_back(new BackwardTr(*b1));
			eraseList.push_back(i);
		}

		/* b->getS1() */
		if (b1->getE1() < b.getS2()) {
			break;
		}
	}

	bool isMatchFound = (overlapList.size() > 0) ? true : false;

	if (isMatchFound) {
		BackwardTr m = b;
		int overlapCount = overlapList.size();
		for (int g = 0; g < overlapCount; g++) {
			// @@ Check
		//	m.merge(overlapList.at(g), m);
			bList->erase(bList->begin() + eraseList.at(g));
		}
		matchBackwardTr(m);
	} else {
		addToBList(new BackwardTr(b));
	
	}

	Util::deleteInVector(&overlapList);
	overlapList.clear();

	eraseList.clear();
}

int DetectorTr::checkMatch(int start, int end) {
	// The vector other is now allocated on the stack
	vector<int> other = vector<int>();
	int lower = scoreList->at(start);
	int upper = scoreList->at(end);
	int otherEnd = end;

	for (int i = start; i <= end; i++) {
		int s = scoreList->at(i);
		if (s != initValue) {
			if (s < lower || s > upper) {
				other.push_back(i);
			}
		}
	}

	int size = other.size();
	if (size > 0) {
		while (size > 0) {
			// The vector coor is now allocated on the stack
			vector<int> coor = vector<int>();
			findTrHelper(other.at(0), &coor);

			int segStart = coor.at(0);
			int segEnd = coor.at(1);
			int segStart2 = scoreList->at(segStart);
			int segEnd2 = scoreList->at(segEnd);
			findMatch(segStart, segEnd, segStart2, segEnd2);

			int segLower = scoreList->at(segStart);
			int segUpper = scoreList->at(segEnd);

			if (segEnd > otherEnd) {
				otherEnd = segEnd;
			}

			// The vector temp is now allocated on the stack
			vector<int> temp = vector<int>();
			for (int j = 0; j < size; j++) {
				int s = scoreList->at(other.at(j));

				if (s != initValue) {
					if (s < segLower || s > segUpper) {
						temp.push_back(other.at(j));
					}
				}
			}
			// Just added start
			other.clear();
			// Just added end
			other = temp;
			size = other.size();
		}
	}

	return otherEnd;
}

vector<BackwardTr *> * DetectorTr::getBList() {
	return bList;
}

} /* namespace tr */
