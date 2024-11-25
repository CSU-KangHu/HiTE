/*
 * FilterTr.cpp
 * 
 *  Created on: Dec 14, 2012
 *      Author: Hani Zakaria Girgis, PhD and Joseph Valencia
 */


#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>

#include "FilterTr.h"
#include "TrKVisitor.h"
#include "TrCsVisitor.h"
#include "TrSineVisitor.h"
#include "../utility/Util.h"
#include "../utility/TSD.h"
#include "../utility/EmptyTSD.h"
#include "../utility/TailFinder.h"
#include "../utility/Tail.h"
#include "../utility/EmptyTail.h"
#include "../utility/GlobAlignE.h"
#include "../utility/LocAlign.h"
#include "../exception/InvalidStateException.h"
#include "../exception/InvalidInputException.h"

using namespace std;
using namespace utility;
using namespace exception;

namespace tr {

FilterTr::FilterTr(string nameIn,const string* seqIn, vector<BackwardTr*>* bListIn, int kIn, int ltrIdIn, int minIn, int maxIn, int minLtrLenIn,int maxLtrLenIn) {
	seq = seqIn;
	cSeq = seq->c_str();
	bList = bListIn;
	name = nameIn;
	k = kIn;
	teList = new vector<LtrTe *>();
	canUseLtr = false;
	canUseLength = true;
	canUseSine = false;
	canUseDNA = true;
	canUsePpt = true;
	canUseTsd = true;

	init = (int) 'N';
	ltrId = ltrIdIn;
   
	tsdW = 20;
	tsdT = 4;

	min = minIn;
	max = maxIn;
	minLtrLen = minLtrLenIn;
	maxLtrLen = maxLtrLenIn;

    tightenBounds();
	//adjust();
	filter();
	removeOverlaps();

	// Sort detections according to the start site
	if(teList->size() > 1){
		sort(teList->begin(), teList->end(), LtrTe::lessThan);
	}
}


FilterTr::~FilterTr() {
	Util::deleteInVector(teList);	
	teList->clear();
	delete teList;
}

void FilterTr::removeOverlaps(){


	if(teList->size() > 0){

		LtrTe * curr = teList->at(0);				

		for(int i = 1; i<teList->size(); i++){
			int currStart = curr->getStart();
			int currEnd = curr->getEnd();

			LtrTe * next = teList->at(i);
			int nextStart = next->getStart();
			int nextEnd = next->getEnd();

			if (currEnd >= nextStart && currEnd <= nextEnd){

				curr->setDeleted(true);
			}
			else{
				curr = next;
			}
		}

		teList->erase(std::remove_if(teList->begin(), teList->end(), [](LtrTe * x){return x->getDeleted();}), teList->end());

	}
}

void FilterTr::tightenBounds()
{   

	int chromLen = seq->length();

	vector<BackwardTr*> * newList = new vector<BackwardTr*>();

	int total = bList->size();
	int kept = 0;

	for (int i=0;i<bList->size();i++)
	{
		BackwardTr * ltr = bList->at(i);

		int s1 = ltr->getS1();
		int e1 = ltr->getE1();
		int upLen = e1+k-s1;

		
		int s2 = ltr->getS2();
		int e2 = ltr->getE2();
		int downLen = e2+k-s2;

		int midpoint = e1+((s2-e1)/2);

 
		float factor = 1.5;

		int smallerLen = std::min(upLen,downLen);
		int largerLen = std::max(upLen,downLen);

		int largerDiff = ceil(factor*largerLen);
		int largerTotal = largerLen+2*largerDiff;

		int smallerDiff = (largerTotal - smallerLen)/2;
		
		int upDiff;
		int downDiff;

		if(smallerLen == upLen){

			upDiff = smallerDiff;
			downDiff = largerDiff;
		}

		else if(smallerLen = downLen){
			downDiff = smallerDiff;
			upDiff = smallerDiff;
		}

		int upTotal = upLen+2*upDiff;
		int downTotal = downLen+2*downDiff;
	
		upDiff = upTotal > maxLtrLen ? (maxLtrLen - upLen)/2 : upDiff;
		downDiff = downTotal > maxLtrLen ? (maxLtrLen-downLen)/2 : downDiff;

		int upStart = s1 - upDiff;
		int upEnd = upStart + upLen + 2 *upDiff + k;

		int downStart = s2 - downDiff;
		int downEnd = downStart + downLen + 2*downDiff + k;

		if (upStart >=0 && upEnd <chromLen && downStart >= 0 && downEnd<chromLen){
			upLen = upEnd > midpoint ? midpoint - upStart : upEnd-upStart; // Does alignment window encroach on midpoint between LTRs?


			string upstream = seq->substr(upStart, upLen);


			downStart = downStart < midpoint ? midpoint : downStart;



			string downstream = seq->substr(downStart, downEnd-downStart);


	
			LocAlign *local = new LocAlign(upstream.c_str(), 0, upstream.length(), downstream.c_str(), 0, downstream.length(), 2, -3, 5, 2);
	
			int newS1 = upStart + local->getQueryStart();
			int newE1 = upStart + local->getQueryEnd();

			int newS2 = downStart + local->getReferenceStart();
			int newE2 = downStart + local->getReferenceEnd();

			double id = local->getIdentity() * 100;

			BackwardTr *updated = new BackwardTr(newS1, newE1, newS2, newE2);

			if (id >= ltrId)
			{  
				updated->setIdentity(id);
				kept++;
				newList->push_back(updated);
			}

			else{
				//cout<<"Removing "<<updated->toString()<<" with %id: "<<id<<endl;
			}



			bool overlaps = true;

			int j = i + 1;

			// Avoiding duplicates
			while (overlaps && j < bList->size())
			{

				BackwardTr *next = bList->at(j);

				int nextStart = next->getS1();

				if (nextStart > updated->getE1())
				{
					overlaps = false;
				}
				else    
				{
					j += 1;
				}
			}
		}
		

	}

	//Removing duplicates
    std::sort( newList->begin(),newList->end(),[](BackwardTr * first, BackwardTr * second){ return first->getS1() < second->getS1();});
	Util::deleteInVector(bList);

	for( int i = 0; i<(int)newList->size()-1;i++){

		BackwardTr * first = newList->at(i);
		int j = i+1;

		BackwardTr * second = newList->at(j);


			while (second->getS1() < first->getE1() && j<newList->size()-1)
		{

			j+=1;
			second = newList->at(j);

		}
		
		bList->push_back(first);
		bList->push_back(second);
		i = j;
	}

	delete newList;
	
}



void FilterTr::adjust() {
	int seqEnd = seq->size() - 1;
	TrKVisitor * kVisitor = new TrKVisitor(k, seqEnd); 
	int size = bList->size();
	for (int i = 0; i < size; i++) {
		bList->at(i)->accept(kVisitor); //adds k to index of each end TR
	}
}


void FilterTr::filter() {
	if (canUseLtr) {
		filterAcc2Ltr();
	} else {
		fillTeList();
	}
	if (canUseSine) {
		filterAcc2Sine();
	}
	
	if(canUseLength)
	{
		filterAcc2Length();
	}

	if (canUseDNA)
	{
		filterAcc2DNA();
	}

	if (canUsePpt)
	{
		filterAcc2Ppt();
	}

	if (canUseTsd) {
		filterAcc2Tsd();
	}
}

void FilterTr::fillTeList() {

	if (teList->size() != 0)
	{ 
			string msg("The TE list must be empty. The current size is: ");
		msg.append(Util::int2string(teList->size()));
		msg.append(".");
		throw InvalidStateException(msg);
	}

	int size = bList->size();

	for (int i = 0; i < size; i++) {
		teList->push_back(
				new LtrTe(bList->at(i), EmptyTSD::getInstance(),
						EmptyTail::getInstance()));
	}
}

void FilterTr::filterAcc2Ltr() {
	if (teList->size() != 0) {
		string msg("The TE list mBackwardTr * ltr = bList->at(i);ust be empty. The current size is: ");
		msg.append(Util::int2string(teList->size()));
		msg.append(".");
		throw InvalidStateException(msg);
	}
	int size = bList->size();

	for (int i = 0; i < size; i++) {
		BackwardTr * ltr = bList->at(i);

		// Overlap between the two LTRs
		int s1 = ltr->getS1();
		int e1 = ltr->getE1();
		int s2 = ltr->getS2();
		int e2 = ltr->getE2();

		if (!Util::isOverlapping(s1, e1, s2, e2)) {

			
			// Distance between the two LTRs
			int sep = abs(s2 - e2 + 1);
			if (sep >= ltrSep) {
				// Length and identity between the two LTRs
				TrCsVisitor * v = new TrCsVisitor(cSeq, minLtrLen, ltrId);
				ltr->accept(v);

				if (v->getIsGood()) {
					teList->push_back(
							new LtrTe(ltr, EmptyTSD::getInstance(),
									EmptyTail::getInstance()));

				}
				else{

				}
				
				delete v;
			}
		}
		
	}
}

void FilterTr::filterAcc2Length(){

	vector<LtrTe *> *temp = new vector<LtrTe *>();
	int size = teList->size();
	int total =size;
	int kept = 0;

	for (int i = 0; i < size; i++)
	{
		LtrTe *te = teList->at(i);

		BackwardTr *ltr = te->getLtr();
		
		int upLen = ltr->getE1()- ltr->getS1();
		int downLen = ltr->getE2() - ltr->getS2();

		int total = ltr->getE2()-ltr->getS1();
		
		bool upFits = upLen >= minLtrLen && upLen <= maxLtrLen;
		bool downFits = downLen >= minLtrLen && downLen <= maxLtrLen;
		bool totalFits = total >= min && total <= max;


		if(upFits && downFits && totalFits){
			LtrTe * longEnough = new LtrTe(*te);
			kept++;
			temp->push_back(longEnough);
		}
		else{

		}
	}

	Util::deleteInVector(teList);
	teList->clear();
	teList = temp;
}

void FilterTr::filterAcc2Sine() {
	vector<LtrTe *> * temp = new vector<LtrTe *>();
	int size = teList->size();

	for (int i = 0; i < size; i++) {
		LtrTe * te = teList->at(i);
		
		BackwardTr * ltr = te->getLtr();
		TrSineVisitor * sineV = new TrSineVisitor(seq, 8, 100, tsdT);
		ltr->accept(sineV);
		bool isTwoSines = sineV->isTwoSines();
		delete sineV;

		if (isTwoSines == false) {
			LtrTe * teWoSine = new LtrTe(*te);
			temp->push_back(teWoSine);
		}
	}

	Util::deleteInVector(teList);
	teList->clear();
	teList = temp;
}


void FilterTr::filterAcc2Tsd() {
	vector<LtrTe *> * temp = new vector<LtrTe *>();

	int size = teList->size();
	for (int i = 0; i < size; i++) {
		LtrTe * te = teList->at(i);	
		BackwardTr * ltr = te->getLtr();
		TSD * tsd = new TSD(seq, ltr, tsdW, init);

		if (tsd->getTsdSize() > tsdT) {
			LtrTe * teWithTsd = new LtrTe(*te);
			teWithTsd->setTsd(tsd);
			temp->push_back(teWithTsd);		
		}
		else{
			te->setTsd(EmptyTSD::getInstance() );
			temp->push_back(new LtrTe(*te));
		}
		delete tsd;
	}

	Util::deleteInVector(teList);
	teList->clear();
	teList = temp;
}

// calculate search window based on minimum and size of interior
int FilterTr::calculateTailWindow(double ratio, int lengthElement, int minimum){	
	int limit = lengthElement > minimum ? minimum : lengthElement;
	int scaled = ceil(ratio*lengthElement);
	return scaled > limit ? scaled : limit;
}

void FilterTr::filterAcc2Ppt() {
	vector<LtrTe *> * temp = new vector<LtrTe *>();
	int size = teList->size();

	for (int i = 0; i < size; i++) {

		LtrTe * te = teList->at(i);
		BackwardTr * ltr = te->getLtr();

		// int minInterior = 50; /*HZG commented out this line*/
		int minWindow = 100; /*HZG moved this line from the if statement below and changed the condition*/
		int interiorLen = ltr->getS2()-ltr->getE1();
		// HZG: Think about this condition
		if(interiorLen >= minWindow){

			Location * location = new Location(ltr->getE1()+1,ltr->getS2()-1);

			int interiorLen = location->getLength();
			double ratio = 0.05;

			int window = calculateTailWindow(ratio,interiorLen,minWindow);

			int seedLen = 5;
			int gapLen = 3;

			TailFinder * finder = new TailFinder(seq, location, TailFinder::MARK_P,seedLen,gapLen, window, tailT);		
			vector<int> * tailVec = finder->getTail();

			if (tailVec->size() == 4) {
				string strand;
				int strandInt = tailVec->at(3);
				if (strandInt == 1) {
					strand = "+";
				} else if (strandInt == -1) {
					strand = "-";
				} else {
					string msg("Invalid strand. ");
					msg.append("Valid strands are 1 and -1, but received: ");
					msg.append(Util::int2string(strandInt));
					msg.append(".");
					cerr <<"About to throw!!!!"<<endl;
					throw InvalidInputException(msg);
				}
				Tail * tail = new Tail(tailVec->at(0), tailVec->at(1), strand, tailVec->at(2));
				LtrTe * teWithTail = new LtrTe(*te);
				teWithTail->setPpt(tail);
				
				delete tail;
				temp->push_back(teWithTail);

			}
			else{
				te->setPpt( EmptyTail::getInstance());
				temp->push_back(new LtrTe(*te));
			}

			delete finder;
			delete location;
		}else{
			te->setPpt( EmptyTail::getInstance());
			temp->push_back(new LtrTe(*te));
		}
	}

	Util::deleteInVector(teList);
	teList->clear();
	teList = temp;
}

void FilterTr::filterAcc2DNA(){

	vector<LtrTe *> *filtered = new vector<LtrTe *>();

	int size = teList->size();
	
	int total = size;
	int kept = 0;

	int minDNASize = 1000;

	for (int i = 0; i < size; i++)
	{

		LtrTe *te = teList->at(i);
		BackwardTr *ltr = te->getLtr();

		//Is there a TIR in the left detection

		int s1 = ltr->getS1();
		int e1 = ltr->getE1();

		int len1 = e1-s1;

		int s2 = ltr->getS2();
		int e2 = ltr->getE2();

		int len2 = e2-s2;



		if(len2<minDNASize && len1 <minDNASize){

		int TIRwindow =30 ;

		len1 = len1< 2*TIRwindow ? len1/2: TIRwindow;

		string deleteThis = seq->substr(s1,len1);



		const char * leftFirst = deleteThis.c_str();

      
		string * temp = new string();
		
		
		Util::revCompDig(cSeq,e1-len1,e1,temp);

		const char * rightFirstRC = temp->c_str();

		LocAlign *firstAlign = new LocAlign(leftFirst, 0, deleteThis.length()- 1, rightFirstRC, 0, temp->length() - 1, 2, -3, 5, 2);

		int minAlignLen = 15;

		float minId = 0.80;

		float id1 = firstAlign->getIdentity();
		

		bool leftTIR = firstAlign->getLength() >= minAlignLen;

		len2 = len2 < 2*TIRwindow ? len2/2 : TIRwindow;

		string deleteThis2 = seq->substr(s2,len2);

		const char *leftSecond = deleteThis2.c_str();

		string * temp2 = new string();

		Util::revCompDig(cSeq, e2 - len2, e2, temp2);

		const char *rightSecondRC = temp2->c_str();

		LocAlign * secondAlign = new LocAlign(leftSecond, 0, deleteThis2.length()-1, rightSecondRC, 0, temp2->length()-1, 2, -3, 5, 2);

		float id2 = secondAlign->getIdentity();

		bool rightTIR = secondAlign->getLength() >= minAlignLen;

		bool twoDNATransposons = leftTIR && rightTIR;

        if(!twoDNATransposons){
			LtrTe * notDNA = new LtrTe(ltr,EmptyTSD::getInstance(),EmptyTail::getInstance());
			kept++;
			filtered->push_back(notDNA);
		}
		else{

		}
		}

		else{
			LtrTe * notDNA = new LtrTe(ltr,EmptyTSD::getInstance(),EmptyTail::getInstance());
			filtered->push_back(notDNA);
		}
	}

		Util::deleteInVector(teList);
		teList->clear();
		teList = filtered;
}

vector<LtrTe *> * FilterTr::getTeList() {
	return teList;
}



string FilterTr::convertNucleotides(string str){
    string ans ="";
	for (int i = 0;i<str.length();i++){
        
		int curr = (int)str.at(i);
        int bp;

		switch(curr){

			case 0:
			    bp = 'A';
				break;
			case 1:
			    bp = 'C';
				break;
			case 2:
				bp = 'G';
				break;
			case 3:
				bp = 'T';
				break;
			default:
			    bp = 'N';

		}

		ans+=bp;

	}

    return ans;

}
void FilterTr::fullFormat(int start, int end)
{   
	ofstream output;
	output.open(bedFileName);

	int step =20;

		for (int i = 0; i < teList->size(); i++)
	{
		LtrTe *curr = teList->at(i);
		BackwardTr *ltr = curr->getLtr();
		
		int s1 = ltr->getS1();
		int e1 = ltr->getE1();
		int len1 = e1-s1;
		int x = maxLtrLen-len1;

		int s2 = ltr->getS2();
		int e2 = ltr->getE2();
		int len2 = e2-s2;
		int y = maxLtrLen-len2;

		int bound = x<y ? x:y;

	    if(len2<minLtrLen && bound<0.75*maxLtrLen){

		GlobAlignE *align = new GlobAlignE(cSeq, s1, e1, cSeq, s2, e2, 1, -1, 4, 1);

		int id = 100* align->getIdentity();

        //expand downstream
		for (int i =0;i<bound;i+=step){

            GlobAlignE * extension = new GlobAlignE(cSeq,e1,e1+step-1,cSeq,e2,e2+step-1,1,-1,5,2);
			int extensionId = extension->getIdentity() *100;

			if(abs(id-extensionId)>=10){
				e2+=step-1;
				e1+=step-1;

				break;
			}

			else if (abs(id - extensionId) >= 10)
			{
	
				break;
			}
		}
        //expand upstream
		for (int i = 0; i < bound; i += step)
		{
            // cout<<"expanding upstream"<<endl;
			GlobAlignE *extension = new GlobAlignE(cSeq, s1-step, s1, cSeq, s2-step, s2, 2, -3, 4, 1);
			int extensionId = extension->getIdentity() * 100;

			if (abs(id - extensionId) >= 10 && extensionId>50)
			{
				s2 -=  step - 1;
				s1 -=  step - 1;
				break;
			}

			else if (abs(id - extensionId) >= 10)
			{
				break;
			}
			
		}

		
		}

		if (e2 - s1 > minLtrLen)
		{
			// HZG: Does this condition need to be here?
		}
	}
	output.close();
}
} /* namespace tr */
