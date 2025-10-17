/*
 * TSD.h
 *
 *  Created on: Dec 21, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef TSD_H_
#define TSD_H_

#include "ITSD.h"
#include "ILocation.h"

namespace utility {

class TSD: public ITSD {
private:
	ILocation * ltTsd;
	ILocation * rtTsd;
	int tsdSize;

public:
	TSD(const string *, ILocation *, int, int);
	TSD(ITSD&);
	TSD(ITSD&,int);
	virtual ~TSD();

	virtual ILocation * getLtTsd();
	virtual ILocation * getRtTsd();
	virtual int getTsdSize();
	virtual string toString();
};

} /* namespace utility */
#endif /* TSD_H_ */
