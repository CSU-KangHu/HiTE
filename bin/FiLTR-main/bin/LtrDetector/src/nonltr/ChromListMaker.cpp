/*
 * ChromListMaker.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: Hani Zakaira Girgis
 *      Modified by Alfredo Velasco II
 */

#include "ChromListMaker.h"

namespace nonltr {

ChromListMaker::ChromListMaker(string seqFileIn) {
	seqFile = seqFileIn;
	chromList = new vector<Chromosome *>();
	chromOList = new vector<ChromosomeOneDigit *>();
	chromSplitMap = new unordered_map<Chromosome *, pair<string, int>>();
	chromOSplitMap =
			new unordered_map<ChromosomeOneDigit *, pair<string, int>>();
	limit = 0;
	passedLimit = false;
}

ChromListMaker::ChromListMaker(string seqFileIn, int limitIn) {
	seqFile = seqFileIn;
	chromList = new vector<Chromosome *>();
	chromOList = new vector<ChromosomeOneDigit *>();
	chromSplitMap = new unordered_map<Chromosome *, pair<string, int>>();
	chromOSplitMap =
			new unordered_map<ChromosomeOneDigit *, pair<string, int>>();

	limit = limitIn;
	if (limit < 1) {
		cerr << "Limit must be positive!" << endl;
		cerr << "`" << limit << "`" << " is invalid." << endl;
		throw std::exception();
	}
	passedLimit = false;
}

ChromListMaker::~ChromListMaker() {
	Util::deleteInVector(chromList);
	delete chromList;
	Util::deleteInVector(chromOList);
	delete chromOList;
	delete chromSplitMap;
	delete chromOSplitMap;
}

const vector<Chromosome *> * ChromListMaker::makeChromList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;
	bool passedLimit = false;
	Chromosome * chrom;
	int seqLen = 0;
	int seqIndex = 0;
	while (in.good()) {
		string line;
		getline(in, line);

		if (line[0] == '>') {
			if (!isFirst) {
				if (passedLimit) {
					passedLimit = false;
				}
				pair<string, int> splitRegion(chrom->getHeader(), seqIndex);
				seqLen = 0;
				seqIndex = 0;
				chrom->finalize();
				if (chrom->getEffectiveSize() > 0) {
					chromList->push_back(chrom);
					chromSplitMap->emplace(chrom, splitRegion);
				} else {
					delete chrom;
				}
			} else {
				isFirst = false;
			}

			chrom = new Chromosome();
			chrom->setHeader(line);
		} else {
			seqLen += line.size();

			if (limit && seqLen > limit) {
				passedLimit = true;
				while (limit && seqLen > limit) {
					const string * base = chrom->getBase();
					string lineSubStr = line.substr(0,
							line.size() - seqLen + limit);
					line = line.substr(line.size() - seqLen + limit);

					chrom->appendToSequence(lineSubStr);
					string oldHeader = chrom->getHeader();
					int baseSize = chrom->getBase()->size();

					pair<string, int> splitRegion(chrom->getHeader(), seqIndex);

					chrom->finalize();

					if (chrom->getEffectiveSize() > 0
							|| chrom->getSegment()->size() != 0) {
						chromList->push_back(chrom);
						chromSplitMap->emplace(chrom, splitRegion);
					} else {
						delete chrom;
					}
					chrom = new Chromosome();
					chrom->setHeader(oldHeader);
					seqLen = line.size();
					seqIndex += baseSize;
				}
				if (line.size() > 0) {
					chrom->appendToSequence(line);
				}

			} else {
				chrom->appendToSequence(line);
			}
		}
		line.clear();
	}
	if (passedLimit) {
		passedLimit = false;
	}
	pair<string, int> splitRegion(chrom->getHeader(), seqIndex);
	chrom->finalize();

	if (chrom->getEffectiveSize() > 0 || chrom->getSegment()->size() != 0) {
		chromList->push_back(chrom);
		chromSplitMap->emplace(chrom, splitRegion);
	} else {
		delete chrom;
	}

	in.close();

	return chromList;
}


const vector<ChromosomeOneDigit *> * ChromListMaker::makeChromOneDigitList() {
	ifstream in(seqFile.c_str());
	bool isFirst = true;

	ChromosomeOneDigit * chrom;
	int seqLen = 0;
	int seqIndex = 0;
	while (in.good()) {
		string line;
		getline(in, line);
		if (line[0] == '>') {
			if (!isFirst) {
				if (passedLimit) {
					passedLimit = false;
				}

				pair<string, int> splitRegion(chrom->getHeader(), seqIndex);
				seqLen = 0;
				seqIndex = 0;
				chrom->finalize();
				if (chrom->getEffectiveSize() > 0) {
					chromOList->push_back(chrom);
					chromOSplitMap->emplace(chrom, splitRegion);
				} else {
					delete chrom;
				}
			} else {
				isFirst = false;
			}

			chrom = new ChromosomeOneDigit();
			chrom->setHeader(line);
		} else {
			seqLen += line.size();
			if (limit && seqLen > limit) {

				passedLimit = true;
				while (limit && seqLen > limit) {
					const string * base = chrom->getBase();
					string lineSubStr = line.substr(0,
							line.size() - seqLen + limit);
					line = line.substr(line.size() - seqLen + limit);

					chrom->appendToSequence(lineSubStr);

					string oldHeader = chrom->getHeader();

					pair<string, int> splitRegion(chrom->getHeader(), seqIndex);
					chrom->finalize();
					int baseSize = chrom->getBase()->size();

					if (chrom->getEffectiveSize() > 0
							|| chrom->getSegment()->size() != 0) {
						chromOList->push_back(chrom);
						chromOSplitMap->emplace(chrom, splitRegion);
					} else {
						delete chrom;
					}

					chrom = new ChromosomeOneDigit();
					chrom->setHeader(oldHeader);
					seqLen = line.size();
					seqIndex += baseSize;
				}

				if (line.size() > 0) {
					chrom->appendToSequence(line);
				}

			} else {
				chrom->appendToSequence(line);
			}
		}
	}
	if (passedLimit) {
		passedLimit = false;
	}
	pair<string, int> splitRegion(chrom->getHeader(), seqIndex);
	chrom->finalize();

	if (chrom->getEffectiveSize() > 0 || chrom->getSegment()->size() != 0) {
		chromOList->push_back(chrom);
		chromOSplitMap->emplace(chrom, splitRegion);
	} else {
		delete chrom;
	}

	in.close();
	return chromOList;
}

pair<string, int> ChromListMaker::getStartOfChromosome(Chromosome * chrom) {
	if (chromSplitMap->find(chrom) != chromSplitMap->end()) {
		return chromSplitMap->at(chrom);
	} else {
		cerr << "Chromosome is not in the list!" << endl;
		throw std::exception();
	}
}
pair<string, int> ChromListMaker::getStartOfChromosome(
		ChromosomeOneDigit * chrom) {
	if (chromOSplitMap->find(chrom) != chromOSplitMap->end()) {
		return chromOSplitMap->at(chrom);
	} else {
		cerr << "Chromosome is not in the list!" << endl;
		throw std::exception();
	}
}

}
/* namespace nonltr */
