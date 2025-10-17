/*
 * BackwardTr.h
 *
 *  Created on: Dec 12, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef BACKWARDTR_H_
#define BACKWARDTR_H_

#include "Tr.h"
#include "ForwardTr.h"
#include "../exception/InvalidStateException.h"
#include "../exception/InvalidOperationException.h"
#include "../utility/Util.h"

namespace tr {
class ForwardTr;

class BackwardTr: public tr::Tr {

public:
	BackwardTr();
	BackwardTr(int, int, int, int);
	virtual ~BackwardTr();
	BackwardTr(BackwardTr&);


	void flip(ForwardTr&);

	bool isOverlapping(BackwardTr *);
	bool isOverlapping(ForwardTr *);

	virtual void initialize(int, int, int, int);

	virtual void checkState();

	virtual int getEnd() const;
	virtual int getStart() const;
	virtual void setEnd(int);
	virtual void setStart(int);
	virtual int getLength();
	virtual string toString();
};
}

#endif /* BACKWARDTR_H_ */
