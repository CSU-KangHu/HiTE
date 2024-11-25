/*
 * ITail.h
 *
 *  Created on: Dec 31, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef ITAIL_H_
#define ITAIL_H_

#include <string>
#include "ILocation.h"

using namespace std;

namespace utility {

class ITail : public utility :: ILocation {
public:
	
	//ITail();
	//virtual ~ITail();
	// Inherited from ILocation
	virtual int getEnd() const = 0;
	virtual int getStart() const= 0;
	virtual void setEnd(int)= 0;
	virtual void setStart(int)= 0;
	virtual int getLength()= 0;
	virtual string toString()= 0;

	// Methods specific to tail objects.
	virtual double getPercentage() const= 0;
	virtual void setPercentage(double)= 0;
	virtual string getStrand() const= 0;
	virtual void setStrand(string)= 0;

};

} /* namespace utility */
#endif /* ITAIL_H_ */
