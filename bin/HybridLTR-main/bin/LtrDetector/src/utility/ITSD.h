/*
 * ITSD.h
 *
 *  Created on: Dec 31, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef ITSD_H_
#define ITSD_H_

#include <string>
#include "ILocation.h"

using namespace std;

namespace utility {

class ITSD {
public:
	virtual ILocation * getLtTsd() = 0;
	virtual ILocation * getRtTsd() = 0;
	virtual int getTsdSize() = 0;
	virtual string toString() = 0;
};

} /* namespace utility */
#endif /* ITSD_H_ */
