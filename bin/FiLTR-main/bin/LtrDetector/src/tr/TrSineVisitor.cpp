/*
 * TrSineVisitor.cpp
 *
 *  Created on: Feb 7, 2013
 *      Author: Hani Zakaria Girgis, PhD
 */
/*
 The purpose of this visitor is to distinguish between true LTR's and two consecutive Sines.
 The existence of TSD's along either side of the two peaks would indicate that they are Sines
 rather than LTRs
*/
#include <iostream>

#include "TrSineVisitor.h"
#include "../utility/Location.h"
#include "../utility/TailFinder.h"
#include "../utility/TSD.h"
#include <algorithm>
#include <cmath>

using namespace utility;
using namespace std;

namespace tr {

TrSineVisitor::TrSineVisitor(const string * seqIn, int tailTIn, int tsdWIn,
		int tsdTIn) {
	seq = seqIn;
	tailT = tailTIn;
	tsdW = tsdWIn;
	tsdT = tsdTIn;

	foundTwoSines = false;
}

TrSineVisitor::~TrSineVisitor() {

}

/**
 * ToDo: check if the length of the LTR is similar to that of SINE < 500bp
 */
void TrSineVisitor::visit(Tr* ltr) {

	int s1 = ltr->getS1();
	int e1 = ltr->getE1();
	int s2 = ltr->getS2();
	int e2 = ltr->getE2();

	// Make locations
	Location * loc1 = new Location(s1, e1);
	Location * loc2 = new Location(s2, e2);

	// Find PloyA tail in the first LTR
	//int win1 = (e1 - s1 + 1) / 2;  // checking the first half 

	int win1 = calculateTailWindow(0.02,loc1->getLength(),50);

	cout<<"PolyA window 1: "<<win1<<endl;

	int seedLen = 5;
	int gapLen = 2;


	TailFinder * f1 = new TailFinder(seq, loc1, TailFinder::MARK_A, seedLen,gapLen, win1,
			tailT);

	// Find PolyA tail in the second LTR
	int win2 = calculateTailWindow(0.02,loc2->getLength(),50);
	cout<<"PolyA window 2: "<<win2<<endl;
	
	TailFinder * f2 = new TailFinder(seq, loc2, TailFinder::MARK_A,seedLen,gapLen, win2,
			tailT);

	// Make sure that the two tails are on the same strand
	if (f1->isTailFound() && f2->isTailFound()) {

		vector <int> * first = f1->getTail();
		vector <int> * second = f2->getTail();

		cout<<"Tail1 ="<< first->at(0)<<":"<<first->at(1)<<endl;;
		cout<<"Tail2 ="<< second->at(0)<<":"<<first->at(1)<<endl;
		/*
		if (f1->getTail()->at(3) == f2->getTail()->at(3)) {
			// Find the first TSD
			TSD * t1 = new TSD(seq, loc1, tsdW, (int) 'N');
			TSD * t2 = new TSD(seq, loc2, tsdW, (int) 'N');
			foundTwoSines = (t1->getTsdSize() > tsdT)
					&& (t2->getTsdSize() > tsdT);
			delete t1;
			delete t2;
		}*/
		foundTwoSines = true;
	}

	delete f1;
	delete f2;

	delete loc1;
	delete loc2;
}
//returns true if this sequence contains two sines rather than an ltr
bool TrSineVisitor::isTwoSines() {
	return foundTwoSines;
}


// calculate search window based on minimum and size of interior
int TrSineVisitor::calculateTailWindow(double ratio, int lengthElement, int minimum){
	// cout << "length interior: " << lengthElement << endl;	
	int limit = lengthElement > minimum ? minimum : lengthElement;
	int scaled = ceil(ratio*lengthElement);
	return scaled > limit ? scaled : limit;
}


/*
 bool TrSineVisitor::isLtrSine(int s, int e) {
 // Make a location object
 Location * l = new Location(s, e);

 // Search for Poly-A tail
 int tailW = (e - s + 1) / 2;
 TailFinder * f = new TailFinder(seq, l, TailFinder::MARK_A, tailW, tailT);
 bool isTailFound = f->isTailFound();

 // Search for TSD
 TSD * t = new TSD(seq, l, tsdW, (int) 'N');
 bool isTsdFound = t->getTsdSize() > tsdT;

 // Free memory
 delete t;
 delete f;
 delete l;
	if (f1->isTailFound() && f2->isTailFound()) {

 // Combine results
 bool result = isTailFound && isTsdFound;
 return result;
 }
 */

} /* namespace tr */
