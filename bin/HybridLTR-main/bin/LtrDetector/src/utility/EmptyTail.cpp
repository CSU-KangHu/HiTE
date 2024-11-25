/*
 * EmptyTail.cpp
 *
 *  Created on: Dec 28, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "EmptyTail.h"
#include "../exception/InvalidOperationException.h"

using namespace exception;

namespace utility {

EmptyTail * EmptyTail::INSTANCE = new EmptyTail();

EmptyTail::EmptyTail() {
	msg = new string("Empty Tail does not allow this operation.");
}

EmptyTail::~EmptyTail() {
	delete msg;
}

EmptyTail* EmptyTail::getInstance() {
	return INSTANCE;
}

string EmptyTail::toString() {
	return string("Empty");
}

int EmptyTail::getEnd() const {
	throw InvalidOperationException(*msg);
}

int EmptyTail::getStart() const {
	throw InvalidOperationException(*msg);
}

void EmptyTail::setEnd(int int1) {
	throw InvalidOperationException(*msg);
}

void EmptyTail::setStart(int int1) {
	throw InvalidOperationException(*msg);
}

int EmptyTail::getLength() {
	throw InvalidOperationException(*msg);
}

double EmptyTail::getPercentage() const {
	throw InvalidOperationException(*msg);
}

void EmptyTail::setPercentage(double double1) {
	throw InvalidOperationException(*msg);
}

string EmptyTail::getStrand() const {
	throw InvalidOperationException(*msg);
}

void EmptyTail::setStrand(string allocator) {
	throw InvalidOperationException(*msg);
}

} /* namespace tr */
