/*
 * Tr.h
 *
 *  Created on: Dec 12, 2012
 *      Author: Hani Zakaria Girgis. PhD
 */

#ifndef TR_H_
#define TR_H_

#include <string>
#include "../utility/ILocation.h"

using namespace std;
using namespace utility;

namespace tr {

class ITrVisitor;

class Tr: public ILocation {

protected:
	int s1;
	int e1;
	int s2;
	int e2;
	bool isOverlappingSame(Tr *);
	bool isOverlappingOpposite(Tr *);
	double id;


public:
	Tr();
	Tr(int, int, int, int);
	Tr(Tr&);
	virtual ~Tr();
	int getS1();
	int getE1();
	int getS2();
	int getE2();


	void setS1(int);
	void setE1(int);
	void setS2(int);
	void setE2(int);

	void accept(ITrVisitor *);
	virtual void initialize(int, int, int, int);

	virtual void checkState() = 0;

	string toString() = 0;
	virtual int getEnd() const = 0;
	virtual int getStart() const = 0;
	virtual void setEnd(int) = 0;
	virtual void setStart(int) = 0;
	virtual int getLength() = 0;
	virtual int getIdentity();
	virtual void setIdentity(int);
};

}

#endif /* TR_H_ */
