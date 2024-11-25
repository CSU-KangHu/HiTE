/*
 * BackwardTr.cpp
 *
 *  Created on: Dec 12, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "BackwardTr.h"


using namespace exception;	
using namespace utility;

namespace tr {

BackwardTr::BackwardTr() :
		Tr() {
}

BackwardTr::BackwardTr(int s1In, int e1In, int s2In, int e2In) :
		Tr(s1In, e1In, s2In, e2In) {
	checkState();
}

BackwardTr::~BackwardTr() {
}

BackwardTr::BackwardTr(BackwardTr& o) :
		Tr(o) {
	checkState();
}

void BackwardTr::checkState() {
	if (s2 < s1) {
		string msg("Invalid direction: s2 must be greater than s1.");
		msg.append(" s1: ");
		msg.append(Util::int2string(s1));
		msg.append(" s2: ");
		msg.append(Util::int2string(s2));
		throw InvalidStateException(msg);
	}
}

void BackwardTr::initialize(int s1, int e1, int s2, int e2) {
	Tr::initialize(s1, e1, s2, e2);
	checkState();
}

/**
 * r: is the flipped TR
 */
void BackwardTr::flip(ForwardTr& r) {
	r.initialize(s2, e2, s1, e1);
}

bool BackwardTr::isOverlapping(BackwardTr* b) {
	return isOverlappingSame(b);
}

bool BackwardTr::isOverlapping(ForwardTr* f) {
	return isOverlappingOpposite(f);
}

int BackwardTr::getStart() const {
	return s1;
}

int BackwardTr::getEnd() const {
	return e2;
}

void BackwardTr::setStart(int sIn) {
	string msg("Setting the start of a TR instance is not allowed.");
	throw InvalidOperationException(msg);
}

void BackwardTr::setEnd(int eIn) {
	string msg("Setting the end of a TR instance is not allowed.");
	throw InvalidOperationException(msg);
}

int BackwardTr::getLength() {
	return e1 - s2 + 1;
}

string BackwardTr::toString() {
	string msg("LTRs_Boundaries ");
	msg.append(Util::int2string(getStart()));
	msg.append("-");
	msg.append(Util::int2string(getEnd()));

	msg.append(" L_LTR ");
	msg.append(Util::int2string(s1));
	msg.append("-");
	msg.append(Util::int2string(e1));

	// msg.append(" L_LTR_Len ");
	// msg.append(Util::int2string((*e2) - (*s2) + 1));

	msg.append(" R_LTR ");
	msg.append(Util::int2string(s2));
	msg.append("-");
	msg.append(Util::int2string(e2));

	// msg.append(" R_LTR_Len ");
	// msg.append(Util::int2string((*e1) - (*s1) + 1));

	return msg;
}
int Tr::getIdentity(){
	return id;
}

void Tr::setIdentity(int newId){
	id = newId;
}

}
