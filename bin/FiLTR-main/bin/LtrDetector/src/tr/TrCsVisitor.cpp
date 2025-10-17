/*
 * TrCsVisitor.cpp
 *
 *  Created on: Dec 14, 2012
 *      Author: Hani Zakaria Girgis, PhD
 */
//Finds Longest common subsequence
#include "TrCsVisitor.h"
#include "Tr.h"
#include "../utility/LCS.h"
#include "../utility/GlobAlignE.h"

// ToDo: delete iostream and print statements
#include <iostream>

using namespace utility;

namespace tr {

TrCsVisitor::TrCsVisitor(const char * seqIn, int minLenIn, int minIdIn) {
	seq = seqIn;
	isGood = false;
	minLen = minLenIn;
	minId = minIdIn;
}

TrCsVisitor::~TrCsVisitor() {
	// TODO Auto-generated destructor stub
}

void TrCsVisitor::visit(Tr* tr) {
	int s1 = tr->getS1();
	int e1 = tr->getE1();
	int s2 = tr->getS2();
	int e2 = tr->getE2();

	//LCS * lcs = new LCS(seq, s1, e1, seq, s2, e2);
	//lcs->printLcs();
	//int lcsScore = lcs->getLenCS();

	double l1 = e1 - s1 + 1;
	//double id1 = 100.00 * (double) lcsScore / l1;

	double l2 = e2 - s2 + 1;
	//double id2 = 100.00 * (double) lcsScore / l2;
    //cout<<"Entering alignment"<<endl;
	GlobAlignE * align = new GlobAlignE(seq, s1,e1,seq,s2,e2,2,-3,5,2);

	double id = 100* align->getIdentity();
	cout<<"ID="<<id<<endl;
	//cout<<"Identity:"<<id<<endl;
	//cout <<align->getIdentity()<<endl;
	//cout<<"l1= "<<l1<<" l2 = "<<l2<<endl;
	if( l1>=minLen && l2 >=minLen && id >=minId){
		isGood = true;
	}

	/*if ((l1 >= minLen && id1 >= minId) || (l2 >= minLen && id2 >= minId)) {
		isGood = true;

		// ToDo: delete print statements
		// Testing start
		/*
		cerr << "L1: " << l1 << " L2: " << l2 << " ID1: " << id1 << " ID2: "
				<< id2 << endl;
		for(int i = s1; i < e1; i++){
			cout << (int) seq[i];
		}
		cout << endl;

		for(int i = s2; i < e2; i++){
			cout << (int) seq[i];
		}
		cout << endl;
		lcs->printLcs();
		cout << "= = =" << endl;
		
		// Testing end
		
	}
*/
	//delete lcs;
	delete align;
}

bool TrCsVisitor::getIsGood() {
	return isGood;
}

} /* namespace tr */
