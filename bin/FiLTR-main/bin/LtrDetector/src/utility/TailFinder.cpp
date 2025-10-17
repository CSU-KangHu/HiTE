/*
 * Assembler.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "TailFinder.h"

#include "../exception/InvalidInputException.h"
#include "../exception/InvalidStateException.h"
#include "../utility/Util.h"

#include <vector>

// #include <iostream>

using namespace std;
using namespace exception;

namespace utility {
TailFinder::TailFinder(const string * seqIn, ILocation * locIn, int whichTailIn, int seedSize,int gapSize,
		int winIn, int minLenIn) {
	seq = seqIn;
	loc = locIn;
	whichTail = whichTailIn;
    // This window is at the end of the region of interest

	seedLen = seedSize;
	gapLen = gapSize;
	win = winIn;
	minLen = minLenIn;

	findMark();
}


TailFinder::~TailFinder() {
	tail->clear();
	delete tail;
}

void TailFinder::findMark() {

	int start = loc->getStart();	
	int end   = loc->getEnd();
	int len   = loc->getLength();

	string * detection = new string(seq->begin() + start, seq->begin()+end+1);

	//cout<<"detection->size() ="<<detection->size()<<endl;

	//int searchEnd = detection->size() / 2;
	//searchEnd = (searchEnd > win) ? win : searchEnd;

	// int searchEnd = win;

	vector<int> * pstvTail = new vector<int>(); //stores info on positive tail
    // cout<<"Repeat:" << endl;
	/*// Test start
	for (int h = 0; h < detection->size(); h++)
	{
		cout << static_cast<int>(detection->at(h)) << " ";
	}
	cout << endl;
	// Test end*/
    // cout<<"Tail Details:"<<endl;

	// Local coordinates
	int localStart = len - win;
	int localEnd   = len - 1;

	if (whichTail == MARK_A) {
		// cout<<"Entering MarkA"<<endl;
		findMarkA(detection, pstvTail, localStart, localEnd);
	} else if (whichTail == MARK_P) {
		findMarkP(detection, pstvTail, localStart, localEnd);
	} else {
		string msg("Invalid mark. Valid marks are: ");
		msg += Util::int2string(MARK_A);
		msg += string(" and ");
		msg += Util::int2string(MARK_P);
		msg += string(". Received: ");
		msg += Util::int2string(whichTail);
		msg += ".";
		throw InvalidInputException(msg);
	}

	// Make global coordinates
	//(start, end, ratio)
	if (pstvTail->size() == 3)
	{
		(*pstvTail)[0] = start + pstvTail->at(0);
		(*pstvTail)[1] = start + pstvTail->at(1);
	}

	// Process the negative tail
	vector<int> *ngtvTail = new vector<int>(); //stores info on reverse complement tail
	string * rcDetection = new string();
	Util::revCompDig(detection, rcDetection);
	delete detection;

	if (whichTail == MARK_A)
	{    
		findMarkA(rcDetection, ngtvTail, localStart,localEnd);
	}
	else if (whichTail == MARK_P)
	{
		findMarkP(rcDetection, ngtvTail, localStart,localEnd);
	}
	delete rcDetection;

	if (ngtvTail->size() == 3)
	{
		// Calculate the global coordinates with respect to the positive strand
		int ee = start + len - ngtvTail->at(0) - 1; //changed from at(1) to reflect alteration of BackwardTr
		int ss = start+len-ngtvTail->at(1) -1;
		//cout << "ss=" << ss+1 << endl;
		//cout << "ee=" << ee << endl;

		(*ngtvTail)[0] = ss;
		(*ngtvTail)[1] = ee;
	}

	// Free memory used by tail(s)
	// The reverse sign is correct
	if (pstvTail->size() == 3 && ngtvTail->size() == 0)
	{
		tail = pstvTail;
		// This -1 means the positive strand?
		tail->push_back(1);
		delete ngtvTail;
	}
	else if (pstvTail->size() == 0 && ngtvTail->size() == 3)
	{
		tail = ngtvTail;
		tail->push_back(-1);
		delete pstvTail;
	}
	else if (pstvTail->size() == 0 && ngtvTail->size() == 0)
	{
		tail = pstvTail;
		delete ngtvTail;
	}
	else if (pstvTail->size() == 3 && ngtvTail->size() == 3)
	{ 
		//if both exist, use the longer one
		int pstvLen = pstvTail->at(1) - pstvTail->at(0) + 1;
		int ngtvLen = ngtvTail->at(1) - ngtvTail->at(0) + 1;
		if (pstvLen > ngtvLen)
		{
			tail = pstvTail;
			tail->push_back(1);
			delete ngtvTail;
		}
		else
		{
			tail = ngtvTail;
			tail->push_back(-1);
			delete pstvTail;
		}
	}
	else
	{
		string msg = string("The tail must have three coordinates only. ");
		msg += string("The +ve tail has ");
		msg += Util::int2string(pstvTail->size());
		msg += string(" coordinates. The -ve tail has ");
		msg += Util::int2string(ngtvTail->size());
		msg += string(" coordinates.");
		throw InvalidStateException(msg);
	}

			// For testing		int winIn, int minLenIn) {

			/*
	 cout << start << "-" << end << endl;

	 if (tail->size() == 4) {
	 cout << ">> Tail length is: " << tail->at(1) - tail->at(0) + 1
	 << " start: " << tail->at(0) << " end: " << tail->at(1)
	 << " percentage is: " << tail->at(2) << "% strand is: "
	 << tail->at(3) << endl;
	 for (int i = start; i <= end; i++) {
	 if (i == tail->at(0)) {
	 cout << "_";
	 }

	 cout << (int) seq->at(i);
	 if (i == tail->at(1)) {
	 cout << "_";		int winIn, int minLenIn) {

	 }
	 }
	 cout << endl;

	 } else {
	 cout << "No PPT was found" << endl;
	 for (int i = start; i <= end; i++) {
	 cout << (int) seq->at(i);
	 }
	 cout << endl;
	 }
	 cout << endl;
	 */
			// End testing
}


string utility::TailFinder::prettyFormatChrom( string * detection){
    std:: string ans;
	for( int i =0;i<detection->length();i++){
		ans+= (int)detection->at(i);
	}
	return ans;
}


//used to find poly(A)
void TailFinder::findMarkA(string * detection, vector<int> * tail, int searchStart,
		int searchEnd) {
	//int code = 3;//code for T nucleotide, which is complement of A in mRNA

	int code = 0; // code for A
	double ratio = 0.70;

	//cout<<"seedLen: "<<seedLen<<endl;
	//cout<<"gapLen: "<<gapLen<<endl;

	// Find the first seed available then immediately break
	bool isSeedFound = false;                                                                
	int seedEnd;


	for (int y = searchStart; y <= searchEnd - seedLen + 1; y++) { 
		int count = 0;
		for (int x = y; x <= y + seedLen - 1; x++) {
			if (detection->at(x) == code) {
				count++;
			}
		}

		if (count == seedLen) {
			seedEnd = y + seedLen - 1;
			isSeedFound = true;
			//cout<<"Seed found at"<<seedEnd-seedLen +1<<","<<seedEnd<<endl;
			break;
		}
	}
    
	// Extend initial seed, allowing for gaps of up to gapLen
	if (isSeedFound) {
		int gapCounter = 0;
		int i = seedEnd + 1;
		for (; i < searchEnd; i++) {
			if (static_cast<int>(detection->at(i)) == code) {
				gapCounter = 0;
			} else {
				gapCounter++;
			}

			if (gapCounter == gapLen) {
				break;
			}
		}
		i--; //decrement to make last valid index before gap
        
		int tailStart = seedEnd - seedLen + 1;
		int tailEnd = i;
		//cout<<"start:"<<tailStart<<" end:"<<tailEnd<<endl;
		//retreats back to last A nucleotide
		for (int j = i; j >= i - gapLen; j--) {
			if (detection->at(j) == code) {
				tailEnd = j;
				break;
			}
		}
        //must end with A
		if (detection->at(tailEnd) != code) {
			string msg("Invalid tail. The tail does not end with T.");
			throw InvalidStateException(msg);
		}
        //cout<<"tailStart:"<<tailStart<<" tailEnd:"<<tailEnd<<endl;
		
	// Extend seed backward

		if(tailStart>searchStart){

		int x = tailStart -1;
		gapCounter = 0;


		for(x; x >searchStart;x--){

			if (detection->at(x) == code) {
					gapCounter = 0;
				} else {
					gapCounter++;
				}
				if (gapCounter == gapLen) {
					break;
				}
		}
//cout << 3 << endl;


		//trims the beginning

		for (int q = x; q <x +gapLen;q++){
			

			if(detection->at(q) == code){

				tailStart = q;
				break;
			}
		}

	}

	
		//count ratio and length

				
		double aTailLen = tailEnd - tailStart + 1;
		double aCount = 0;

		//counts A
		for (int h = tailStart; h <= tailEnd; h++) {
			if (detection->at(h) == code) {
				aCount++;
			}
		}

		if (aTailLen >= minLen && (aCount / aTailLen) >= ratio) {
			tail->push_back(tailStart); //tail[0] = start
			tail->push_back(tailEnd); //tail[1] = end
			//double tCount = 0;
			//cout<<"tCount1="<<tCount<<endl;
	

			tail->push_back(100 * aCount / aTailLen);
		}
/*
		if (tail->size() == 0 && searchEnd - searchStart + 1 >= minLen) { //there is enough left to fit in a minLen tail
			int shift = searchStart + seedEnd; // increment by 
			string * rest = new string(detection->begin() + seedEnd,
					detection->begin() + detection->size());
			findMarkA(rest, tail, shift, searchEnd - seedEnd);
			delete rest;
		}*/
	}
}

//used to find purines (A and G)
void TailFinder::findMarkP(string * detection, vector<int> * tail, int searchStart,
		int searchEnd) {

	double ratio = 0.70;
	int codeA = 0;
	int codeG = 2;

	// Find a seed
	bool isSeedFound = false;
	int seedEnd;


	//cout<<"Starting purine search at: "<<searchStart<<endl;

	
	for (int y = searchStart; y <= searchEnd - seedLen + 1; y++) {
		int count = 0;
		for (int x = y; x <= y + seedLen - 1; x++) {
			if (detection->at(x) == codeG || detection->at(x) == codeA) {
				count++;
			}
		}

		if (count == seedLen) {
			seedEnd = y + seedLen - 1;
			isSeedFound = true;
			break;
		}
	}

	// Extend 
	if (isSeedFound) {
		int gapCounter = 0;
		int i = seedEnd + 1;

		//Extend forward
		for (; i < searchEnd; i++) {
			if (detection->at(i) == codeG || detection->at(i) == codeA) {
				gapCounter = 0;
			} else {
				gapCounter++;
			}
			if (gapCounter == gapLen) {
				break;
			}
		}

		i--; //decrement to make last valid index before gap

//cout << 1 << endl;

		// Trims end of the tail
		int tailStart = seedEnd - seedLen + 1;
		int tailEnd = i;

		//cout<<"TAIL_START "<<tailStart<<endl;
		//cout<<"TAIL_END "<<tailEnd<<endl;

		for (int j = i; j >= i - gapLen; j--) {
			if (detection->at(j) == codeG || detection->at(j) == codeA) {
				tailEnd = j;
				break;
			}
		}

//cout << 2 << endl;

		// Post condition to ensure that the tail ends with G or A
		if (!(detection->at(tailEnd) == codeG || detection->at(tailEnd) == codeA)) {
			string msg("The tail does not end with C or T. ");
			msg.append("The invalid base is: ");
			msg.append(Util::int2string((int) detection->at(tailEnd)));
			msg.append(".");
			throw InvalidStateException(msg);
		}

		// Extend seed backward

		if(tailStart>searchStart){

		int x = tailStart -1;
		gapCounter = 0;

		

		for(x; x >searchStart;x--){

			if (detection->at(x) == codeG || detection->at(x) == codeA) {
					gapCounter = 0;
				} else {
					gapCounter++;
				}
				if (gapCounter == gapLen) {
					break;
				}
		}
//cout << 3 << endl;


		//trims the beginning

		for (int q = x; q <x +gapLen;q++){
			

			if(detection->at(q) == codeG || detection->at(q) == codeA){

				tailStart = q;
				break;
			}
		}

	}

//cout << 4 << endl;
		// Count AG in the tail
		// cout << "End of tail: " << tailEnd << endl;
		double pTailLen = tailEnd - tailStart + 1;
		double agCount = 0;
		for (int h = tailStart; h <= tailEnd; h++) {
			if (detection->at(h) == codeG || detection->at(h) == codeA) {
				agCount++;
			}
		}

		// Verifies the ratio and length criteria
		if (pTailLen >= minLen && (agCount / pTailLen) >= ratio) {
			tail->push_back(tailStart );
			tail->push_back(tailEnd);

			/*double pCount = 0;
			for (int g = tailStart; g <= tailEnd; g++) {
				if (detection->at(g) == codeG || detection->at(g) == codeA) {
					pCount++;
				}
			}*/
			tail->push_back(100 * agCount / pTailLen);
		}

		double count = 100*agCount/pTailLen;

		//cout<<"Found percentage "<<count<<endl;
		//cout<<"Found length "<<pTailLen<<endl;

		
		// recursively find next seed after the first one
		if (tail->size() == 0 && searchEnd - searchStart + 1 >= minLen) {
			// int shift = seedEnd;
			// string * rest = new string(detection->begin() + seedEnd,
			//			detection->begin() + detection->size());
			//cout<<"seedEnd: "<<seedEnd<<endl;
			//cout<<"searchEnd: "<<searchEnd<<endl;
			findMarkP(detection, tail, seedEnd, searchEnd);
		}
	}
}

/*
 * The size of the tail indicates the following:
 * 0: no tail is found
 * 4: start, end, strand: 1 indicates pstv and -1 indicates ngtv.
 */
vector<int> * TailFinder::getTail() {

	if(tail->size()==4){
		/*
		cout<< "Chromosome length is: "<<seq->length()<<endl;
		for(int x = 0;x<tail->size();x++){
			cout<<tail->at(x)<<",";

		}
		cout<<endl;
		*/

		int start = tail->at(0);
		int end = tail->at(1);
		int strand = tail->at(3);

		/*

		if(strand ==-1 || end<start){

			int temp = start;
			start = end;
			end = temp;
			cout<<"Swapped 'em"<<endl;
		}

		cout<<"s: "<<start<<endl;
		cout<<"e: "<<end<<endl;*/

		
		string * detection = new string(seq->begin() + start, seq->begin() + end + 1);

			for (int h = 0; h < detection->size(); h++){
				
			//	cout << static_cast<int>(detection->at(h));
			}
		//	cout << endl;

	
	}
	return tail;
}

/**
 * If the vector representing the tail has a size of zero,
 * no tail is found.
 */
bool TailFinder::isTailFound() {
	bool r = tail->size() == 0 ? false : true;
	return r;
}

}
