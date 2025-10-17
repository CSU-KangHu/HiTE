/*
 * TrPptVisitor.cpp
 *
 *  Created on: Dec 26, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "TrPptVisitor.h"
#include "../exception/InvalidStateException.h"

namespace tr {

TrPptVisitor::TrPptVisitor(const string * seqIn) {
	seq = seqIn;
	win = 500;
}

TrPptVisitor::~TrPptVisitor() {

}

void TrPptVisitor::visit(Tr * trIn) {
	tr = trIn;
	segStart = tr->getStart();
}

void TrPptVisitor::	searchPstv(){
	int (* cTable) = new int[win][2];


}

} /* namespace tr */
