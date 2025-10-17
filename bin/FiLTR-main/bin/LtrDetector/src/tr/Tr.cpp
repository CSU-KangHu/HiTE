/*
 * Tr.cpp
 *
 *  Created on: Dec 12, 2012
 *      Author: girgishz
 */

#include "Tr.h"
#include "ITrVisitor.h"
#include "../exception/InvalidStateException.h"
#include "../utility/Util.h"

// #include <stdlib.h>
// ToDo: delete
#include <iostream>
// end delete

using namespace std;
using namespace utility;
using namespace exception;

namespace tr {

/**
 * No-parameter constructor
 */
Tr::Tr() {
	s1 = 0;
	e1 = 0;
	s2 = 0;
	e2 = 0;
}

Tr::Tr(int s1In, int e1In, int s2In, int e2In) {
	initialize(s1In, e1In, s2In, e2In);
}

Tr::Tr(Tr& o) {
	initialize(o.getS1(), o.getE1(), o.getS2(), o.getE2());

	id = o.getIdentity();
}

void Tr::initialize(int s1In, int e1In, int s2In, int e2In) {
	s1 = s1In;
	e1 = e1In;
	s2 = s2In;
	e2 = e2In;
	id = 0;

	if (e1In < s1In) {
		string msg("The start of the first TR must be <= its end. ");
		msg.append("The start is: ");
		msg.append(Util::int2string(s1In));
		msg.append(" The end is: ");
		msg.append(Util::int2string(e1In));
		throw InvalidStateException(msg);
	}

	if (e2In < s2In) {
		string msg("The start of the second TR must be <= its end. ");
		msg.append("The start is: ");
		msg.append(Util::int2string(s2In));
		msg.append(" The end is: ");
		msg.append(Util::int2string(e2In));
		throw InvalidStateException(msg);
	}

	// checkState();
}

Tr::~Tr() {
}

int Tr::getS1() {
	return s1;
}

int Tr::getS2() {
	return s2;
}

int Tr::getE1() {
	return e1;
}

int Tr::getE2() {
	return e2;
}

void Tr::setS1(int s1In) {
	s1 = s1In;
}

void Tr::setS2(int s2In) {
	s2 = s2In;
}

void Tr::setE1(int e1In) {
	e1 = e1In;
}

void Tr::setE2(int e2In) {
	e2 = e2In;
}

bool Tr::isOverlappingSame(Tr* same) {
	bool cond1 = Util::isOverlapping(s1, e1, same->getS1(), same->getE1());
	bool cond2 = Util::isOverlapping(s2, e2, same->getS2(), same->getE2());
	return cond1 && cond2;
}

bool Tr::isOverlappingOpposite(Tr* opposite) {
	bool cond1 = Util::isOverlapping(s1, e1, opposite->getS2(),
			opposite->getE2());
	bool cond2 = Util::isOverlapping(s2, e2, opposite->getS1(),
			opposite->getE1());

	return cond1 && cond2;
}

void Tr::accept(ITrVisitor * v) {
	v->visit(this);
}

}
