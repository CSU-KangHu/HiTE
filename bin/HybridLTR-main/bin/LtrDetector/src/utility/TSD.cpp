/*
 * TSD.cpp
 *
 *  Created on: Dec 21, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "TSD.h"
#include "Location.h"
#include "EmptyLocation.h"
#include "Util.h"
#include "LCSubStr.h"
#include "../exception/InvalidStateException.h"

#include <vector>

using namespace std;
using namespace exception;

namespace utility {

TSD::TSD(const string * seq, ILocation * te, int w, int init) {
	const char * cSeq = seq->c_str();

	// Determine the left window
	int ltEnd = te->getStart() - 1;
	if (ltEnd < 0) {
		ltEnd = 0;
	}
	int ltStart = ltEnd - w + 1;
	if (ltStart < 0) {
		ltStart = 0;
	}
	int temp = ltEnd;
	for (int i = ltEnd; i >= ltStart; i--) {
		if (cSeq[i] == init) {
			break;
		} else {
			temp = i;
		}
	}
	ltStart = temp;
	Location * ltWin = new Location(ltStart, ltEnd);

	// Determine the right window
	int rtStart = te->getEnd() + 1;
	int rtEnd = rtStart + w - 1;
	int lastIndex = seq->size() - 1;
	if (rtEnd > lastIndex) {
		rtEnd = lastIndex;
	}
	int temp1 = rtStart;
	for (int i = rtStart; i <= rtEnd; i++) {
		if (cSeq[i] == init) {
			break;
		} else {
			temp1 = i;
		}
	}
	rtEnd = temp1;
	Location * rtWin = new Location(rtStart, rtEnd);

	// Determine the closest TSD
	LCSubStr * lcss = new LCSubStr(cSeq, ltWin, cSeq, rtWin);

	vector<vector<ILocation *> *> * r = lcss->getCSubStr();

	int min = 1000000;
	int minIndex = -1;

	int tsdFound = r->at(0)->size();

	for (int j = 0; j < tsdFound; j++) {
		
		ILocation * lt = r->at(0)->at(j);
		ILocation * rt = r->at(1)->at(j);

		int ltDis = te->getStart() - lt->getEnd();
		int rtDis = rt->getStart() - te->getEnd();
		int dis = ltDis + rtDis;

		if (ltDis < 0 || rtDis < 0) {
			string msg("Distance cannot be negative. The left distance is: ");
			msg.append(Util::int2string(ltDis));
			msg.append(" The right distance is: ");
			msg.append(Util::int2string(rtDis));
			throw InvalidStateException(msg);
		}

		if (dis < min) {
			min = dis;
			minIndex = j;
		}
	}

	int ltSize;
	int rtSize;
	if (minIndex == -1) {
		ltTsd = EmptyLocation::getInstance();
		rtTsd = EmptyLocation::getInstance();
		ltSize = 0;
		rtSize = 0;

	} else {
		ltTsd = new Location(*(r->at(0)->at(minIndex)));
		rtTsd = new Location(*(r->at(1)->at(minIndex)));
		ltSize = ltTsd->getLength();
		rtSize = rtTsd->getLength();

	}

	if (ltSize != rtSize) {
		string msg("The two sites must have the same length. ");
		msg.append("The length of the left site is: ");
		msg.append(Util::int2string(ltSize));
		msg.append(" The length of the right site is: ");
		msg.append(Util::int2string(rtSize));
		msg.append(".");
		throw InvalidStateException(msg);
	}
	tsdSize = ltSize;

	// Free resources
	delete lcss;
	delete ltWin;
	delete rtWin;
}

TSD::TSD(ITSD& copy) {
	ltTsd = new Location(*copy.getLtTsd());
	rtTsd = new Location(*copy.getRtTsd());
	tsdSize = copy.getTsdSize();
}

TSD::TSD(ITSD& copy, int offset) {


	ltTsd = new Location(copy.getLtTsd()->getStart()+offset,copy.getLtTsd()->getEnd()+offset);

	rtTsd = new Location(copy.getRtTsd()->getStart()+offset,copy.getRtTsd()->getEnd()+offset);
	tsdSize = copy.getTsdSize();
}

TSD::~TSD() {
	if (ltTsd != EmptyLocation::getInstance()) {
		delete ltTsd;
	}

	if (rtTsd != EmptyLocation::getInstance()) {
		delete rtTsd;
	}
}

ILocation * TSD::getLtTsd() {
	return ltTsd;
}

ILocation * TSD::getRtTsd() {
	return rtTsd;
}

int TSD::getTsdSize() {
	return tsdSize;
}

string TSD::toString() {
	string msg("L_TSD ");
	msg.append(ltTsd->toString());
	msg.append(" R_TSD ");
	msg.append(rtTsd->toString());
	return msg;
}

} /* namespace utility */
