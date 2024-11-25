/*
 * InvalidOperationException.h
 *
 *  Created on: Dec 20, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#ifndef INVALIDOPERATIONEXCEPTION_H_
#define INVALIDOPERATIONEXCEPTION_H_

#include <string>
#include <stdexcept>

using namespace std;

namespace exception {

class InvalidOperationException : public std::runtime_error{
public:
	InvalidOperationException(string msg);
	//virtual ~InvalidOperationException();
};

}

#endif /* INVALIDOPERATIONEXCEPTION_H_ */
