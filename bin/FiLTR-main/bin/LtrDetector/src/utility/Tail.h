/*
 * Tail.h
 *
 *  Created on: Dec 27, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef TAIL_H_
#define TAIL_H_

#include <string>
#include "ITail.h"

using namespace std;

namespace utility {

class Tail: public utility::ITail {
private:
	int start;
	int end;
	double percentage;
	string strand;
	void initialize(int, int, string, double);

public:
	Tail(int, int, string, double);
	Tail(ITail&);
	Tail(ITail&,int);
	virtual ~Tail();

	// Inherited from ILocation
	virtual int getEnd() const;
	virtual int getStart() const;
	virtual void setEnd(int);
	virtual void setStart(int);
	virtual int getLength();
	virtual string toString();

	// Methods specific to tail objects.
	virtual double getPercentage() const;
	virtual void   setPercentage(double);
	virtual string getStrand() const;
	virtual void setStrand(string);
};

} /* namespace tr */
#endif /* TAIL_H_ */
