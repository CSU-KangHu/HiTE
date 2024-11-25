/*
 * Tail.cpp
 *
 *  Created on: Dec 27, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include <string>
#include "Tail.h"
#include "Util.h"

namespace utility {

Tail::Tail(int startIn, int endIn, string strandIn, double percentageIn) {
	initialize(startIn, endIn, strandIn, percentageIn);
}

Tail::Tail(ITail& copy) {
	initialize(copy.getStart(), copy.getEnd(), copy.getStrand(),
			copy.getPercentage());
}

Tail::Tail(ITail& copy,int offset) {
	initialize(copy.getStart()+offset, copy.getEnd()+offset, copy.getStrand(),
			copy.getPercentage());
}

void Tail::initialize(int startIn, int endIn, string strandIn,
		double percentageIn) {
	start = startIn;
	end = endIn;
	strand = strandIn;
	percentage = percentageIn;
}

Tail::~Tail() {
	// TODO Auto-generated destructor stub
}

int Tail::getStart() const {
	return start;
}

void Tail::setStart(int startIn) {
	start = startIn;
}

int Tail::getEnd() const {
	return end;
}

void Tail::setEnd(int endIn) {
	end = endIn;
}

int Tail::getLength() {
	return end - start + 1;
}

string Tail::toString() {
	string msg ="";
	msg.append(Util::int2string(start));
	msg.append("\t");
	msg.append(Util::int2string(end));
	msg.append("\t");
	msg.append(strand);
	msg.append("\t");
	msg.append(Util::double2string(percentage));
	return msg;
}

double Tail::getPercentage() const {
	return percentage;
}

void Tail::setPercentage(double percentageIn) {
	percentage = percentageIn;
}

string Tail::getStrand() const {
	return strand;
}

void Tail::setStrand(string strandIn) {
	strand = strandIn;
}

} /* namespace tr */
