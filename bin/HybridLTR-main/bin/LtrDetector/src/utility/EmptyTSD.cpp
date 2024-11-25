/*
 * EmptyTSD.cpp
 *
 *  Created on: Dec 28, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "Util.h"
#include "EmptyTSD.h"
#include "../exception/InvalidOperationException.h"

using namespace exception;

namespace utility {

EmptyTSD * EmptyTSD::INSTANCE = new EmptyTSD();

EmptyTSD::EmptyTSD() {
	msg = new string("Empty TSD does not allow this operation.");
}

EmptyTSD::~EmptyTSD() {
	delete msg;
}

string EmptyTSD::toString() {
	return string("Empty");
}

EmptyTSD * EmptyTSD::getInstance() {
	return INSTANCE;
}

ILocation* EmptyTSD::getLtTsd() {
	throw InvalidOperationException(*msg);
}

ILocation* EmptyTSD::getRtTsd() {
	throw InvalidOperationException(*msg);
}

int EmptyTSD::getTsdSize() {
	throw InvalidOperationException(*msg);
}

} /* namespace utility */
