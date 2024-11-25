/*
 * ITrVisitor.h
 *
 *  Created on: Dec 14, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef ITRVISITOR_H_
#define ITRVISITOR_H_

namespace tr {

class Tr;

class ITrVisitor {
public:
	virtual void visit(Tr *) = 0;
};

} /* namespace tr */
#endif /* ITRVISITOR_H_ */
