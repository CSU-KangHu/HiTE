/*
 * EmptyTSD.h
 *
 *  Created on: Dec 28, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef EMPTYTSD_H_
#define EMPTYTSD_H_

#include "ITSD.h"
#include "ILocation.h"

namespace utility {

class EmptyTSD : public ITSD{
private:
	string * msg;
	static EmptyTSD * INSTANCE;
	EmptyTSD();

public:
	virtual ~EmptyTSD();
	virtual ILocation * getLtTsd();
	virtual ILocation * getRtTsd();
	virtual int getTsdSize();
	virtual string toString();
	static EmptyTSD * getInstance();
};

} /* namespace utility */
#endif /* EMPTYTSD_H_ */
