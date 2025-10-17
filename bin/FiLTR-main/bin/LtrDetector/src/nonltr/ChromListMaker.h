/*
 * ChromListMaker.h
 *
 *  Created on: Mar 13, 2014
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef CHROMLISTMAKER_H_
#define CHROMLISTMAKER_H_

#include <string>
#include <vector>
#include <unordered_map>

#include "Chromosome.h"
#include "ChromosomeOneDigit.h"

#include "../utility/Util.h"

using namespace std;
using namespace utility;

namespace nonltr {

class ChromListMaker {
private:
	vector<Chromosome *> * chromList;
	vector<ChromosomeOneDigit *> * chromOList;
	unordered_map<Chromosome *, pair<string, int>> * chromSplitMap;
	unordered_map<ChromosomeOneDigit *, pair<string, int>> * chromOSplitMap;
	string seqFile;
	int limit;
	bool passedLimit;

public:
	ChromListMaker(string);
	ChromListMaker(string, int);
	virtual ~ChromListMaker();
	const vector<Chromosome *> * makeChromList();
	const vector<ChromosomeOneDigit *> * makeChromOneDigitList();
	pair<string, int> getStartOfChromosome(Chromosome *);
	pair<string, int> getStartOfChromosome(ChromosomeOneDigit *);
};

} /* namespace nonltr */
#endif /* CHROMLISTMAKER_H_ */
