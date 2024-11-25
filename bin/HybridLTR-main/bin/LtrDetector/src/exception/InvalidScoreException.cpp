/*
 * InvalidScoreException.cpp
 *
 *  Created on: Apr 27, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */

#include "InvalidScoreException.h"

#include <string>
#include <iostream>

using namespace std;
namespace exception{

InvalidScoreException::InvalidScoreException(string massage) {
	cerr << "Invalid Score Exception." << endl;
	cerr << massage << endl;
}

InvalidScoreException::~InvalidScoreException() {
	// TODO Auto-generated destructor stub
}
}
