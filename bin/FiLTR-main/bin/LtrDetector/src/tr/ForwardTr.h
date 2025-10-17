/*
 * ForwardTr.h
 *
 *  Created on: Dec 12, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef FORWARDTR_H_
#define FORWARDTR_H_

#include "Tr.h"

namespace tr {

class BackwardTr;

class ForwardTr: public tr::Tr {
public:
	ForwardTr();
	ForwardTr(int, int, int, int);
	virtual ~ForwardTr();
	ForwardTr(ForwardTr&);

	void merge(BackwardTr *, ForwardTr&);
	void merge(ForwardTr *, ForwardTr&);
	void flip(BackwardTr&);

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

#endif /* FORWARDTR_H_ */
