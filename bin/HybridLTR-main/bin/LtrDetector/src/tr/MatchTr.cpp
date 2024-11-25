/*
  Created by Joseph V
  encia 21 February 2018
*/

#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <stack>
#include <cmath>

#include "MatchTr.h"
#include "ForwardTr.h"
#include "Candidate.h"
#include "PairContainer.h"
#include "../utility/Util.h"
#include "../exception/InvalidInputException.h"

using namespace std;

using namespace utility;
using namespace exception;

namespace tr{

MatchTr::MatchTr(vector<int> * scoreListIn,int kIn, int initValueIn,int minIn, int maxIn, int ltrMinIn,int plateauLenIn,  int gapTolIn, int id)   {
    
    k = kIn;
    scoreList = scoreListIn;
  
    initValue = initValueIn;
    

    min = minIn;
	max = maxIn;
	ltrMin = ltrMinIn;
	
    minPlateauLen = plateauLenIn;
	diffThresh = gapTolIn;
	gapTol = gapTolIn;
	identity = id;

    bList = new vector<BackwardTr *>();
    
    cleanAndMerge();


} 

void MatchTr::cleanAndMerge(){
    
	int len = scoreList->size();

  
	for(int i = 0;i<len;i++){
        
	    int curr = scoreList->at(i);
        int next;
		
		if(curr != initValue){
			
			int peakLength = 1;
			
			for(int j = i+1;j<len;j++){  
				next = scoreList->at(j);
				if(next == initValue){
					break;
				}
			    peakLength++;
				
			}
			if(peakLength < minPlateauLen){
				
				std::tuple<char,int, int> discard = make_tuple('D',i, i + peakLength);

				spikes.push_back(discard);

			}
			else{
				std::tuple<char,int, int> keep = make_tuple('K',i, i + peakLength);

                spikes.push_back(keep);
			}
			
			i += peakLength;
		}
	}

	forwardMerge();


	backwardMerge();

	//medianSmooth();
	smoothAndCollect();


}


void MatchTr::forwardMerge(){
	// Thanks to Robert Hubley for finding and fixing a bug related to the following 3 lines.
	if(spikes.size() ==0){
		return;
	}
    
    for (int i = 0; i < spikes.size()-1; i++)
	{
		char curr_type;
		int curr_start;
		int curr_end;
		std::tie(curr_type, curr_start, curr_end) = spikes.at(i);

		int level = findMedian(curr_start, curr_end);

		char next_type;
        int next_start;
		int next_end;
		std::tie(next_type, next_start, next_end) = spikes.at(i + 1);

		int neighborScore = findMedian(next_start,next_end);

		if (next_type == 'K')
		{
			if (abs(neighborScore - level) < diffThresh && (curr_end + gapTol) >= next_start) // spikes are at same level and within distance of gapTol
			{   
				
				for (int j = curr_end; j <= next_start; j++) //replace curr_end with curr_start
				{
					
					(*scoreList)[j] = neighborScore;
				}
				if (curr_type == 'D')
				{
					spikes[i] = make_tuple('K', curr_start,next_start-1); //flip curr to a section to keep

			
				}
			}

					
		}

			else if (curr_type =='K' && next_type == 'D'){
			 
				if (abs(neighborScore - level) < diffThresh && (curr_end + gapTol) >= next_start) // spikes are at same level and within distance of k

				{
					for (int j = curr_end; j <= next_start; j++) //replaced curr_end with curr_start. TODO: Find out why next-start bounds are off at 1824000-1832000
					{
						
						(*scoreList)[j] = level;
					}
					
						spikes[i + 1] = make_tuple('K', curr_start, next_end); //flip next to a section to keep
				}
				
			}
			
	}

}

int MatchTr::findMedian(int start, int end){
	vector<int> section1(scoreList->begin() + start, scoreList->begin() + end);
	if(section1.size() ==0){
		return 0;
	}
	else if(section1.size() ==1){
		return section1.at(0);
	}
	else if (section1.size()%2 ==0){
	std::nth_element(section1.begin(), section1.begin() + section1.size() / 2, section1.end());
	
	return section1[section1.size() / 2];
	}
	else{
		std::nth_element(section1.begin(), section1.begin() + section1.size() / 2, section1.end());

		return section1[(section1.size() / 2)+1];
	}

}

void MatchTr::medianSmooth()
{

	vector<int> *temp = new vector<int>();

	int step = 20;
	int i = 0;

	while (i < scoreList->size() - step)
	{

		int median = findMedian(i, i + step);

		for (int j = i; j < i + step; j++)
		{
			temp->push_back(median);
		}

		i += step;
	}
	int score = findMedian(i, scoreList->size());
	for (int j = i; j < scoreList->size(); j++)
	{
		temp->push_back(score);
	}

	scoreList->clear();
	scoreList = temp;
}

void MatchTr::backwardMerge()
{
	// Thanks to Robert Hubley for finding and fixing a bug related to the following 3 lines.
	if(spikes.size() ==0){
		return;
	}
    
	for (int i = spikes.size()-1; i>=1; i--)
	{
		char curr_type;
		int curr_start;
		int curr_end;
		std::tie(curr_type, curr_start, curr_end) = spikes.at(i);
		int level = findMedian(curr_start, curr_end);
		
		char next_type;
        int next_start;
		int next_end;
		std::tie(next_type, next_start, next_end) = spikes.at(i - 1);

		int neighborScore = findMedian(next_start,next_end);


		if (curr_type == 'K' && next_type == 'D'){

			if (abs(neighborScore - level) < diffThresh && (curr_start - gapTol) <= next_end) // spikes are at same level and within distance of k
			{
				for (int j = curr_start; j > next_start; j--)
				{
					(*scoreList)[j] = level;
				}
			}

		}
	}
}

void MatchTr::smoothAndCollect(){
	//Repeat collection
	int len = scoreList->size();
	for (int i = 0; i < len; i++)
	{

		int curr = scoreList->at(i);
		int next;

		if (curr != initValue)
		{

			int peakLength = 1; 

			for (int j = i + 1; j < len; j++)
			{
				next = scoreList->at(j);
			
				if (next == initValue)
				{
					break;
				}
				peakLength++;
			}
		
			int minSizeKeep = identity * ltrMin / 100;

			if (peakLength >= minSizeKeep) //added this parameter
			{
                int height = findMedian(i,i+peakLength-1);

				Candidate *keep = new Candidate(i, i + peakLength - 1, height);
				plateaus.push_back(keep);

			}
			else
			{

				for (int k = i; k < i + peakLength; k++)
				{
					(*scoreList)[k] = initValue;
				}

			}
			i += peakLength - 1;
		}
	}

       PairContainer * matcher = new PairContainer(min,max,diffThresh);
	
        Candidate * curr;
		Candidate * match;
		BackwardTr * pair;
		for (int i = 0; i < plateaus.size(); i++)
		{
			curr = plateaus.at(i);
            match = matcher->hashOrReturn(curr);
			if(match!=nullptr){
				pair = new BackwardTr(match->getStart(), match->getEnd(), curr->getStart(), curr->getEnd());
				bList->push_back(pair);

			}
	}
		
} 

bool MatchTr::isMatch(int firstStart,int secondStart){

	int firstHeight = scoreList->at(firstStart);

	int matchLoc = firstStart + firstHeight;

	return matchLoc == secondStart;

}

vector<BackwardTr *> * MatchTr::getRepeatCandidates(){
    return bList;
}

MatchTr::~MatchTr(){
	bList->clear();
	delete bList;
}

vector<int>* MatchTr::getScoreList(){
	return scoreList;
}

}