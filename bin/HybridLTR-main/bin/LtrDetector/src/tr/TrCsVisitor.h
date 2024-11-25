/*
 * TrCsVisitor.h
 *
 *  Created on: Dec 14, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef TRCSVISITOR_H_
#define TRCSVISITOR_H_

#include "ITrVisitor.h"

namespace tr {

class TrCsVisitor: public tr::ITrVisitor {
private:
	const char * seq;
	bool isGood;
	int minLen;
	int minId;

public:
	TrCsVisitor(const char *, int, int);
	virtual ~TrCsVisitor();
	virtual void visit(Tr *);
	bool getIsGood();
};

} /* namespace tr */
#endif /* TRCSVISITOR_H_ */
