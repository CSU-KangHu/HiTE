/*
 * TrCollector.cpp
 *
 *  Created on: Jan 2, 2013
 *      Author: Hani Zakaria Girgis, PhD and Joseph Valencia
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include "TrCollector.h"
#include "ScorerTr.h"
#include "DetectorTr.h"
#include "FilterTr.h"
#include "MatchTr.h"
#include "BackwardTr.h"
#include "LtrTe.h"
#include "../utility/ITail.h"
#include "../utility/ITSD.h"
#include "../utility/EmptyTSD.h"
#include "../utility/EmptyTail.h"
#include "../utility/EmptyLocation.h"
#include "../nonltr/ChromosomeOneDigit.h"

#include "../utility/Util.h"

using namespace std;
using namespace utility;

namespace tr {

TrCollector::TrCollector(ChromosomeOneDigit * chromIn /*HZG changed this parameter*/, std::string bedFileNameIn, std::string nameIn, int minIn,
						 int maxIn, int ltrMinIn, int ltrMaxIn, int idIn, int kIn, int plateauLenIn, int gapTolIn,bool printRawIn, bool printCleanIn,bool bedFormatIn, bool nestedIn)
{	
	k = kIn;
	min = minIn;
	max = maxIn;
	minPlateauLen = plateauLenIn;

	gapTol = gapTolIn;
	name = nameIn;
	ltrMin = ltrMinIn;
	ltrMax = ltrMaxIn;
	bedFileName = bedFileNameIn;
    identity = idIn;
	chrom = chromIn;
	cout << "Processing " << chrom->getHeader() << " ..." << endl;

	bedFormat = bedFormatIn;
	printClean = printCleanIn;
	printRaw = printRawIn;
	displayNested = nestedIn;

	nestedTeList = new vector<LtrTe*>();

	collect();
}

void TrCollector::collect() {
    int maxScore = max - ltrMin <=2000 ? 2000: max-ltrMin;
	int minScore = min - ltrMax <= 2000 ? 2000 : min - ltrMax;

	ScorerTr * scorer = new ScorerTr(chrom, k, minScore,maxScore);
	int length = chrom->getEffectiveSize();
	
	int init = scorer->getInitialScore();

	vector<int> * scoreList = scorer->getScores();
	

	if(printRaw){

		scoresFormat(scoreList,bedFileName,"Raw");
	}

	

	MatchTr * matcher = new MatchTr(scoreList,k, init, min,max,ltrMin, minPlateauLen, gapTol,identity);



	if(printClean){

		scoresFormat(matcher->getScoreList(),bedFileName,"Clean");
	}
    
	cout << "Filtering " << chrom->getHeader() << " ..." << endl;
	FilterTr *filter = new FilterTr(name, chrom->getBase(), matcher->getRepeatCandidates(), k,identity,min,max,ltrMin,ltrMax);

	teList = filter->getTeList();

	// Write regular TE
	string regularFile = bedFileName+"Detector.bed";
	#pragma omp critical
	{
		outputAnnotation(teList,regularFile);
		cout << "Output from: " << name << " found in: " << regularFile<< endl;
	}	

	// Free memory. The scorer requires large memory
	delete matcher;
	delete scorer;
	

	teList = filter->getTeList();

	if(displayNested){
		findNested();

		// Write nested elements
		#pragma omp critical
		{
			string nestedFile = bedFileName+"NestedDetector.bed";
			outputAnnotation(nestedTeList,nestedFile);

			cout << " Nested output from: " << name << " found in: " << nestedFile << endl;
		}
	}

	delete filter;
}

TrCollector::~TrCollector() {
	// HZG added the following code
	Util::deleteInVector(nestedTeList);
	nestedTeList->clear();
	delete nestedTeList;
}

void TrCollector::scoresFormat(vector<int> * scores,string outputDir,string type){


	ofstream output;

	string scoreOutput = outputDir+type+"Scores.csv";
	
	output.open(scoreOutput);

	for (int i = 0;i<scores->size();i++){

		 output << i << "," << scores->at(i)<< endl;
	}
	output.close();
}

void TrCollector::outputAnnotation(vector<LtrTe*> * myTEList, string fileName)
{
	ofstream output;
	// Append to file
	output.open(fileName, std::fstream::out | std::fstream::app);

	if(output.is_open()){

	int size = myTEList->size();

	if(!bedFormat){

	output<<"SeqID"<<"\t"<<"Retrotransposon"<<"\t"<<"Left_LTR"<<"\t"<<"Right_LTR"<<"\t\t"<<"Left_TSD"<<"\t"<<"Right_TSD"<<"\t"<<"Polypurine Tract"<<"\t\t"<<"TG"<<"\t"<<"CA"<<endl;
	output<<"\t"<<"Start"<<"\t"<<"End"<<"\t"<<"Start"<<"\t"<<"End"<<"\t"<<"Start"<<"\t"<<"End"<<"\t"<<"ID"<<"\t"<<"Start"<<"\t"<<"End"<<"\t"<<"Start"<<"\t"<<"End"<<"\t"<<"Start"<<"\t"<<"End"<<"\t"<<"Strand"<<"\t"<<"Purine\%"<<"\t"<<"Start"<<"\t"<<"End"<<endl;
	}
	for (int i = 0; i < size; i++)
	{
		LtrTe *curr = myTEList->at(i);

		BackwardTr* ltr = curr->getLtr();

		if(bedFormat){
			output << chrom->getHeader().substr(1) << "\t" << curr->getStart() << '\t' << curr->getEnd() << "\t" << ltr->getS1()<< '\t' << ltr->getE1() << "\t" << ltr->getS2() << '\t' << ltr->getE2()<< endl;
		}
		else{
			output << chrom->getHeader().substr(1) << "\t";
			int id = ltr->getIdentity();
			output << curr->getStart()<< '\t' << curr->getEnd() << "\t" << ltr->getS1() << '\t' << ltr->getE1()<<"\t"<<ltr->getS2()<< "\t" <<ltr->getE2()<<"\t"<<id<<"\t" ;

			//NOTE: if either member of pair is string::no_pos then print dash


			ITail * ppt = curr->getPpt();	
			ITSD * tsd = curr->getTsd();


			if(tsd != EmptyTSD::getInstance()){

				ILocation * leftTSD = tsd->getLtTsd();
				ILocation * rightTSD = tsd->getRtTsd();

				output<<leftTSD->toString()<<"\t"<<rightTSD->toString()<<"\t";
			}
			else{
				output<<"---"<<"\t"<<"---"<<"\t"<<"---"<<"\t"<<"---""\t";
			}

			if(ppt != EmptyTail::getInstance()){

				output<<ppt->toString()<<"\t";
			}

			else{

				output<<"---"<<"\t"<<"---"<<"\t"<<"0"<<"\t"<<"--"<<"\t";

			}

			ILocation * tgcaMotif = curr->getTgCaMotif(chrom->getBase());

			if(tgcaMotif!=EmptyLocation::getInstance()){
				output<<tgcaMotif->toString()<<endl;;
			}
			else{
				output<<"---"<<"\t"<<"---"<<endl;
			}
		}
	}
	output.close();
	}
	else{
		cerr<<fileName<<" cannot be opened."<<endl;
		throw std::exception();
	}
}

void TrCollector::findNested() {
	cout << "Finding nested LTR-RT in " << chrom->getHeader() << endl;

	const vector<vector<int> *> * segment = chrom->getSegment();

	int candidateCount = teList->size();

	if (candidateCount > 0) {
		int firstCandIndex = 0;
		int lastCandIndex = 0;
		int segmentNumber = segment->size();
		for (int i = 0; i < segmentNumber; i++) {
			vector<int> * s = segment->at(i);
			ILocation * c = teList->at(firstCandIndex);
			// A segment may have no detections
			if (Util::isOverlapping(s->at(0), s->at(1), c->getStart(),
					c->getEnd())) {
				lastCandIndex = findNestedHelper1(s->at(0), s->at(1),
						firstCandIndex);
				findNestedHelper2(s->at(0), s->at(1), firstCandIndex, lastCandIndex);
				firstCandIndex = lastCandIndex + 1;
				if (firstCandIndex >= candidateCount) {
					break;
				}
			}
		}
	}
}


int TrCollector::findNestedHelper1(int segStart, int segEnd, int firstCandIndex) {
	ILocation * cand = teList->at(firstCandIndex);
	if (!Util::isOverlapping(segStart, segEnd, cand->getStart(),
			cand->getEnd())) {
		string msg("The first candidate is not overlapping with the segment. ");
		msg.append("Candidate location is: ");
		msg.append(cand->toString());
		msg.append(" Segment location is: ");
		msg.append(Util::int2string(segStart));
		msg.append("-");
		msg.append(Util::int2string(segEnd));
		throw InvalidInputException(msg);
	}

	int lastCandIndex = -1;
	int candidateNumber = teList->size();
	for (int c = firstCandIndex; c < candidateNumber; c++) {
		ILocation * cand = teList->at(c);
		if (Util::isOverlapping(segStart, segEnd, cand->getStart(),
				cand->getEnd())) {
			lastCandIndex = c;
		} else {
			break;
		}
	}

	if (lastCandIndex < 0) {
		string msg("The index of the last candidate cannot be negative.");
		throw InvalidStateException(msg);
	}

	return lastCandIndex;
}

void TrCollector::findNestedHelper2(int segStart, int segEnd, int firstCandIndex,
		int lastCandIndex){
	
	const string * base = chrom->getBase();

	// Build a chromosome from the expanded region
	for (int i = firstCandIndex; i <= lastCandIndex; i++){
		LtrTe * inner = teList->at(i);
		int innerSize = inner->getLength();
		int expandedStart = inner->getStart() - max/2;
		int expandedEnd = inner->getEnd() + max/2;
		if(expandedStart < segStart){
			expandedStart = segStart;
		}
		if(expandedEnd > segEnd){
			expandedEnd = segEnd;
		}
 
		string converted = Util::oneDigitToNuc(base->substr(expandedStart, expandedEnd - expandedStart + 1));
		string tempName("temp");
		ChromosomeOneDigit * outerChrom = new ChromosomeOneDigit(converted, tempName);

		int minScore = innerSize+1;
		int maxScore = expandedEnd - expandedStart+1-ltrMin;

		if(innerSize > segEnd - segStart + 1){
			cerr << "\tTrCollector::trainHelper2: " << endl;
			cerr << "\tTE size is greater than segment size." << endl;
			cerr << "\tTE: " << inner->toString () << endl;
			cerr << "\tSegment: " << segStart << "--" <<  segEnd << endl;
			cerr << "\tSkipping this element for now." << endl;
			// HZG: Keep this code.
			/*
			cerr << "Base size: " << base->size() << endl;
			string segment = Util::oneDigitToNuc(base->substr(segStart,segEnd - segStart + 1));
			cerr << "String of segment: " << segment << endl << endl;

			string elementOneDigit = base->substr(inner->getStart(), inner->getEnd() - inner->getStart() + 1);
			string element = Util::oneDigitToNuc(elementOneDigit);
			cerr << "String of element: " << element << endl << endl;
			*/
			// Print the above message for now.
			//throw std::exception();
		} else if(minScore > maxScore){
			cerr << "\tTrCollector::trainHelper2: " << endl;
			cerr << "\tThe maximum distance cannot be less than the minimum distance. " << endl;
			cerr << "\tTE: " << inner->toString () << endl;
			cerr << "\tSegment: " << segStart << "--" <<  segEnd << endl;
			cerr << "\tSkipping this element for now." << endl << endl;
		} else{
			// Make a scorer, a matcher, and a filter for the vacinity of each TE
			ScorerTr * scorer = new ScorerTr(outerChrom,k,minScore,maxScore);
			MatchTr * matcher = new MatchTr(scorer->getScores(),k,scorer->getInitialScore(),minScore,maxScore,ltrMin,minPlateauLen,gapTol,identity);
			FilterTr * filter = new FilterTr(name, outerChrom->getBase(), matcher->getRepeatCandidates(), k, identity,minScore,maxScore,ltrMin,ltrMax);
			
			// Collected nested TE (if any)
			vector<LtrTe*> * nested = filter->getTeList();
			for( int x = 0; x < nested->size(); x++){
				LtrTe * ltr = nested->at(x);
				LtrTe * shifted = new LtrTe(*ltr,expandedStart);
				nestedTeList->push_back(shifted);
			}

			// Free memory
			delete scorer;
			delete matcher;
			delete filter;
		}
		delete outerChrom;
	}
}

void TrCollector::printIndex(string outputFile, vector<LtrTe *> * teList) {
	ofstream outIndex;
	outIndex.open(outputFile.c_str(), ios::out /*| ios::app*/);

	// Write the index of the repeat segment [x,y[ "exclusive" with respect with the start (chrK:start-end)
	string header = chrom->getHeader();
	int size = teList->size();
	for (int j = 0; j < size; j++) {
		LtrTe * te = teList->at(j);
		outIndex << te->toString(header) << endl;
	}
	outIndex.close();
}

void TrCollector::printMasked(string outputFile, vector<LtrTe *> * teList) {
	
	/*string baseCopy = string(*(chrom->getBase()));
	int size = teList->size();

	for (int j = 0; j < size; j++) {
		LtrTe * te = teList->at(j);
		int teStart = te->getStart();
		int teEnd = te->init, k, getEnd();

		for (int h = teStart; h <= teEnd; h++) {
			baseCopy[h] = tolower(baseCopy[h]);
		}
	}

	ofstream outMask;
	outMask.open(outputFile.c_str(), ios::out /*| ios::app);
	outMask << chrom->getHeader() << endl;
	int step = 50;
	int len = baseCopy.size();
	for (int i = 0; i < len; i = i + step) {
		int e = (i + step - 1 > len - 1) ? len - 1 : i + step - 1;
		for (int k = i; k <= e; k++) {
			outMask << baseCopy[k];
		}
		outMask << endl;
	}
	outMask.close();*/
}

} /* namespace tr */
