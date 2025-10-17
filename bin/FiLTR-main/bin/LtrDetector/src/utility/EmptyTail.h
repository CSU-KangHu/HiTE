/*
 * EmptyTail.h
 *
 *  Created on: Dec 28, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef EMPTYTAIL_H_
#define EMPTYTAIL_H_

#include <string>
#include "ITail.h"

using namespace std;

namespace utility {

class EmptyTail: public ITail {
private:
	string * msg;
	static EmptyTail * INSTANCE;
	EmptyTail();

public:
	virtual ~EmptyTail();
	static EmptyTail * getInstance();

	// Inherited from ILocation
	virtual int getEnd() const;
	virtual int getStart() const;
	virtual void setEnd(int);
	virtual void setStart(int);
	virtual int getLength();
	virtual string toString();

	// Methods specific to tail objects.
	virtual double getPercentage() const;
	virtual void setPercentage(double);
	virtual string getStrand() const;
	virtual void setStrand(string);
};

} /* namespace tr */
#endif /* EMPTYTAIL_H_ */
