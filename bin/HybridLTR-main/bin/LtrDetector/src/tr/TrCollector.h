/*
 * TrCollector.h
 *
 *  Created on: Jan 2, 2013
 *      Author: Hani Zakaria Girgis, PhD
 */
// 
#ifndef TRCOLLECTOR_H_
#define TRCOLLECTOR_H_

#include "LtrTe.h"
#include "../nonltr/ChromosomeOneDigit.h"

using namespace nonltr;

namespace tr {

class TrCollector {
private:
	ChromosomeOneDigit * chrom;
	int k;
	int min; // minimum separation distance between ltr
	int max; //absolute maximum ^
	int ltrMin;
	int ltrMax;
	//int d; // "delta" incremental separation between ltr on iteration
	int minPlateauLen;
	int diffThresh;
	
	int gapTol;
	int identity;
	std::string csvFileName;
	std::string bedFileName;
	std::string name;
	// vector<LtrTe *> * teList; /*HZG commented out this file*/
	bool bedFormat;
	bool printRaw;
	bool printClean;
	bool displayNested;

	// A list of regular (very likely not nested) TE
	vector<LtrTe*> * teList;
	// A list of nested TE
	vector<LtrTe*> * nestedTeList;

	void collect();
	void findNested();
	int findNestedHelper1(int, int, int);
	void findNestedHelper2(int, int, int, int);
	void scoresFormat(vector<int>*,string,string);

public:
	TrCollector(ChromosomeOneDigit * /*HZG changed this parameter*/,std::string,std::string, int, int,int,int,int,int,int,int,bool,bool,bool,bool);
	virtual ~TrCollector();
	// vector<LtrTe *> * getTeList();
	void printIndex(string, vector<LtrTe *> * );
	void printMasked(string, vector<LtrTe *> *);
	void outputAnnotation(vector<LtrTe *> *,string);

};

} /* namespace tr */
#endif /* TRCOLLECTOR_H_ */
