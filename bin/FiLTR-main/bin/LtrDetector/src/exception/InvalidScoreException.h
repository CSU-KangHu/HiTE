/*
 * InvalidScoreException.h
 *
 *  Created on: Apr 27, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef INVALIDSCOREEXCEPTION_H_
#define INVALIDSCOREEXCEPTION_H_

#include <string>

using namespace std;

namespace exception{
	class InvalidScoreException {
	public:
		InvalidScoreException(string);
		virtual ~InvalidScoreException();
	};
}

#endif /* INVALIDSCOREEXCEPTION_H_ */
