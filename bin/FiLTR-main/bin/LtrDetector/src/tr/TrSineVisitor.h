/*
 * TrSineVisitor.h
 *
 *  Created on: Feb 7, 2013
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef TRSINEVISITOR_H_
#define TRSINEVISITOR_H_

#include "ITrVisitor.h"
#include "Tr.h"

namespace tr {

class TrSineVisitor: public tr::ITrVisitor {

private:
	const string * seq;
	int tsdW;
	int tsdT;
	int tailT;
	bool foundTwoSines;



public:
	TrSineVisitor(const string *, int, int, int);
	virtual ~TrSineVisitor();
	virtual void visit(Tr *);
	bool isTwoSines();
	int calculateTailWindow(double,int,int);
};

} /* namespace tr */
#endif /* TRSINEVISITOR_H_ */
