/*
 * TrPptVisitor.h
 *
 *  Created on: Dec 26, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef TRPPTVISITOR_H_
#define TRPPTVISITOR_H_

#include "ITrVisitor.h"
#include "../utility/ILocation.h"

using namespace utility;

namespace tr {

class TrPptVisitor: public ITrVisitor {
private:
	const string * seq;
	Tr * tr;
	int segStart;
	int win;

	void searchPstv();

public:
	TrPptVisitor(const string *);
	virtual ~TrPptVisitor();
	virtual void visit(Tr *);
};

} /* namespace tr */
#endif /* TRPPTVISITOR_H_ */
