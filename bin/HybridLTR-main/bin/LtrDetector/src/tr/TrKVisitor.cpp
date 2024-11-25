/*
 * TrKVisitor.cpp
 *
 *  Created on: Dec 14, 2012
 *      Author: zakarota
 */

#include "TrKVisitor.h"

namespace tr {

TrKVisitor::TrKVisitor(int kIn, int endIn) {
	k = kIn;
	end = endIn;
}

TrKVisitor::~TrKVisitor() {
	// TODO Auto-generated destructor stub
}

void TrKVisitor::visit(Tr* tr) {
	int e1 = tr->getE1() + k - 1;
	if (e1 > end) {
		e1 = end;
	}
	tr->setE1(e1);

	int e2 = tr->getE2() + k - 1;
	if (e2 > end) {
		e2 = end;
	}
	tr->setE2(e2);
}

} /* namespace tr */
