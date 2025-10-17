/*
 * TrKVisitor.h
 *
 *  Created on: Dec 14, 2012
 *      Author: zakarota
 */

#ifndef TRKVISITOR_H_
#define TRKVISITOR_H_

#include "Tr.h"
#include "ITrVisitor.h"

namespace tr {

class TrKVisitor: public tr::ITrVisitor {
private:
	int k;
	int end;

public:
	TrKVisitor(int, int);
	virtual ~TrKVisitor();
	virtual void visit(Tr *);
};

} /* namespace tr */
#endif /* TRKVISITOR_H_ */
