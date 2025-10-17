/*
 * ForwardTr.cpp
 *
 *  Created on: Dec 12, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */
#include "Tr.h"
#include "ForwardTr.h"
#include "BackwardTr.h"
#include "../exception/InvalidStateException.h"
#include "../exception/InvalidOperationException.h"
#include "../utility/Util.h"

#include <iostream>
using namespace std;

using namespace exception;
using namespace utility;

namespace tr {

ForwardTr::ForwardTr() :
		Tr() {

}

ForwardTr::ForwardTr(int s1In, int e1In, int s2In, int e2In) :
		Tr(s1In, e1In, s2In, e2In) {
	checkState();
}

ForwardTr::~ForwardTr() {
}

ForwardTr::ForwardTr(ForwardTr& o) :
		Tr(o) {
	checkState();
}

void ForwardTr::checkState() {
	if (s1 > s2) {
		string msg("Invalid direction: s2 must be greater than s1.");
		msg.append(" s1: ");
		msg.append(Util::int2string(s1));
		msg.append(" s2: ");
		msg.append(Util::int2string(s2));
		throw InvalidStateException(msg);
	}
}

void ForwardTr::initialize(int s1, int e1, int s2, int e2) {
	Tr::initialize(s1, e1, s2, e2);
	checkState();
}

/**
 * f: is the other TR to be merged with this TR
 * r: is the resulting TR
 */
void ForwardTr::merge(ForwardTr * f, ForwardTr& r) {
	int s1f = f->getS1();
	int e1f = f->getE1();
	int s2f = f->getS2();
	int e2f = f->getE2();

	int s1n = (s1 < s1f) ? s1 : s1f;
	int e1n = (e1 > e1f) ? e1 : e1f;
	int s2n = (s2 < s2f) ? s2 : s2f;
	int e2n = (e2 > e2f) ? e2 : e2f;

	r.initialize(s1n, e1n, s2n, e2n);
}

/**
 * f: is the other TR to be merged with this TR
 * r: is the resulting TR
 */
void ForwardTr::merge(BackwardTr * b, ForwardTr& r) {
	int s1b = b->getS1();
	int e1b = b->getE1();
	int s2b = b->getS2();
	int e2b = b->getE2();

	int s1n = (s1 < s2b) ? s1 : s2b;
	int e1n = (e1 > e2b) ? e1 : e2b;
	int s2n = (s2 < s1b) ? s2 : s1b;
	int e2n = (e2 > e1b) ? e2 : e1b;

	r.initialize(s1n, e1n, s2n, e2n);
}

/**
 * r: is the flipped TR
 */
void ForwardTr::flip(BackwardTr& r) {
	r.initialize(s2, e2, s1, e1);
}

bool ForwardTr::isOverlapping(BackwardTr* b) {
	return isOverlappingOpposite(b);
}

bool ForwardTr::isOverlapping(ForwardTr* f) {
	return isOverlappingSame(f);
}

int ForwardTr::getStart() const {
	return s1;
}

int ForwardTr::getEnd() const {
	return e2;
}

void ForwardTr::setStart(int sIn) {
	string msg("Setting the start of a TR instance is not allowed.");
	throw InvalidOperationException(msg);
}

void ForwardTr::setEnd(int eIn) {
	string msg("Setting the end of a TR instance is not allowed.");
	throw InvalidOperationException(msg);
}

int ForwardTr::getLength() {
	return e2 - s1 + 1;
}

string ForwardTr::toString() {
	string msg("TwoLTRs ");
	msg.append(Util::int2string(getStart()));
	msg.append("-");
	msg.append(Util::int2string(getEnd()));

	msg.append(" L_LTR ");
	msg.append(Util::int2string(s1));
	msg.append("-");
	msg.append(Util::int2string(e1));

	//msg.append(" L_LTR_Len ");
	//msg.append(Util::int2string((*e1) - (*s1) + 1));

	msg.append(" R_LTR ");
	msg.append(Util::int2string(s2));
	msg.append("-");
	msg.append(Util::int2string(e2));

	//msg.append(" R_LTR_Len ");
	//msg.append(Util::int2string((*e2) - (*s2) + 1));

	return msg;
}

}
