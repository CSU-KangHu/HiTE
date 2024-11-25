/*
 * ScorerTr.cpp
 *
 *  Created on: Nov 30, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "ScorerTr.h"
#include "../nonltr/KmerHashTable.h"
#include "ForwardTr.h"
#include "../nonltr/ChromosomeOneDigit.h"
#include "../utility/Util.h"
#include "../exception/InvalidInputException.h"

#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stack>
#include <algorithm>
using namespace std;

using namespace nonltr;
using namespace utility;
using namespace exception;

namespace tr {

const int ScorerTr::INITIAL_VALUE = -10;
const int ScorerTr::INITIAL_SCORE = 0;

ScorerTr::ScorerTr(ChromosomeOneDigit *chromIn, int motifSizeIn, int minIn,int maxIn)
{
	chrom = chromIn;
	k = motifSizeIn;
	min = minIn;
	max = maxIn;

	if (max < min) {
		string msg(
				"The maximum distance cannot be less than the minimum distance. ");
		msg.append("The minimum distance is: ");
		msg.append(Util::int2string(min));
		msg.append(". The maximum distance is: ");
		msg.append(Util::int2string(max));
		msg.append(".");
		throw InvalidInputException(msg);
	}

	kmerTable = new KmerHashTable<int,int>(k, INITIAL_VALUE);
	scoreList = new vector<int>(chrom->getBase()->size(), INITIAL_SCORE);

	scoreNew();
}

ScorerTr::~ScorerTr() {
	delete kmerTable;
	scoreList->clear();
	delete scoreList;
	
}

void ScorerTr::score() {
	
	const vector<vector<int> *> * segment = chrom->getSegment();
	const char * segBases = chrom->getBase()->c_str();

	for (int s = 0; s < segment->size(); s++) {
		
		int start = segment->at(s)->at(0);
        int end = segment->at(s)->at(1);
		
		vector<int> * hashList = new vector<int>();
		kmerTable->hash(segBases, start, end - k + 1, hashList);

		// I commented out the +1. It does not make sense not to score the first word of a segment.
		for (int i = start ; i <= end - k + 1; i++) {
			int keyHash = hashList->at(i - start);
			int lastIndex = kmerTable->valueOf(keyHash);

			if (lastIndex != INITIAL_VALUE) {
				int d1 = abs(i - lastIndex);
				if (d1 >= min && d1 <= max) {
					(*scoreList)[i] = lastIndex;
					int secondLastIndex = scoreList->at(lastIndex);
					if (abs(lastIndex - secondLastIndex) > d1) {
						(*scoreList)[lastIndex] = i;
					}
				}
			}

			kmerTable->insert(keyHash, i);
		}
		hashList->clear();
		delete hashList;

		// Handle last word
		for (int i = end - k + 2; i <= end; i++) {
			(*scoreList)[i] = scoreList->at(i - 1);
		}
	}

	// Test code
	ofstream output;
	string file1 = "scores.txt";
	output.open(file1);
	cout<<"Before: " << endl;
	for(auto e : *scoreList){
        output << e<< " ";
	}
	output<<endl;
	// End test code

	// Generate distance --- positive and negative --- from indexes.
	for(int i = 0;i<scoreList->size();i++){
		int score = scoreList->at(i);
		if(score != INITIAL_SCORE){
			(*scoreList)[i] = score - i;
		}
	}

	// Test code
	cout << "After: " << endl;
	for (auto e : *scoreList){   
		output <<e<<" ";
	}
	output << endl;
	output.close();
	// End test code
}

void ScorerTr::scoreNew() {
	
	const vector<vector<int> *> * segment = chrom->getSegment();
	const char * segBases = chrom->getBase()->c_str();

	for (int s = 0; s < segment->size(); s++) {
		int start = segment->at(s)->at(0);
		int end = segment->at(s)->at(1);

		vector<int> * hashList = new vector<int>();
		kmerTable->hash(segBases, start, end - k + 1, hashList);

		// I commented out the +1. It does not make sense not to score the first word of a segment.
		for (int i = start /*+ 1*/; i <= end - k + 1; i++) {
			int keyHash = hashList->at(i - start);
			int lastIndex = kmerTable->valueOf(keyHash);
			if (lastIndex != INITIAL_VALUE) { 
				
				int d1 = abs(i - lastIndex);
				if (d1 >= min && d1 <= max)
				{
					(*scoreList)[i] = lastIndex -i;
					int scoreAtLastIndex = scoreList->at(lastIndex); 
					if (scoreAtLastIndex == INITIAL_SCORE || d1 < abs(scoreAtLastIndex))
					{
						(*scoreList)[lastIndex] = i - lastIndex;
					}
				}
			}

			kmerTable->insert(keyHash, i);
		}
		hashList->clear();
		delete hashList;

		// Handle last word
		for (int i = end - k + 2; i <= end; i++) {
			(*scoreList)[i] = scoreList->at(i - 1);
		}
	}
}

int ScorerTr::findMedian(int start, int end)
{
	vector<int> section1(scoreList->begin() + start, scoreList->begin() + end);
	if (section1.size() == 0)
	{
		return 0;
	}
	else if (section1.size() == 1)
	{
		return section1.at(0);
	}
	else if (section1.size() % 2 == 0)
	{
		std::nth_element(section1.begin(), section1.begin() + section1.size() / 2, section1.end());

		return section1[section1.size() / 2];
	}
	else
	{
		std::nth_element(section1.begin(), section1.begin() + section1.size() / 2, section1.end());

		return section1[(section1.size() / 2) + 1];
	}
}

void ScorerTr::medianSmooth()
{

	vector<int> *temp = new vector<int>();

	int step = 20;
	int i = 0;

	while (i < scoreList->size() - step)
	{

		int median = findMedian(i, i + step);

		for (int j =i;j < i + step;j++)
		{
			temp->push_back(median);
		}

		i+=step;
	}
    int score = findMedian(i,scoreList->size());
	for(int j=i;j<scoreList->size();j++){
		temp->push_back(score);
	}

	scoreList->clear();
	scoreList = temp;
}

vector<int>* ScorerTr::getScores() {
	return scoreList;
}

int ScorerTr::getInitialScore() {
	return INITIAL_SCORE;
}

}
