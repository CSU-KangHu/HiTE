/*
 * ChromosomeOneDigit.cpp
 *
 *  Created on: Jul 31, 2012
 *      Author: Hani Zakaria Girgis, PhD at the NCB1/NLM/NIH
 * A	A
 * T	T
 * G	G
 * C	C
 * R	G or A
 * Y	T or C
 * M	A or C
 * K	G or T
 * S	G or C
 * W	A or T
 * H	A or C or T
 * B	G or T or C
 * V	G or C or A
 * D	G or T or A
 * N	G or T or A or C
 */
#include <iostream>
#include <map>

#include "Chromosome.h"
#include "ChromosomeOneDigit.h"
#include "../exception/InvalidInputException.h"

using namespace exception;

namespace nonltr {

ChromosomeOneDigit::ChromosomeOneDigit() :
		Chromosome() {
}

ChromosomeOneDigit::ChromosomeOneDigit(string fileName) :
		Chromosome(fileName) {
	help();
}

ChromosomeOneDigit::ChromosomeOneDigit(string fileName, int segmentLength,
		int maxLength) :
		Chromosome(fileName, segmentLength, maxLength) {
	help();
}

ChromosomeOneDigit::ChromosomeOneDigit(string& seq, string& info) :
		Chromosome(seq, info) {
	help();
}

ChromosomeOneDigit::ChromosomeOneDigit(string& seq, string& info, int length) :
		Chromosome(seq, info, length) {
	help();
}

void ChromosomeOneDigit::help() {
	// Build codes
	buildCodes();
	// Modify the sequence in the super class
	encodeNucleotides();
}

void ChromosomeOneDigit::finalize() {
	Chromosome::finalize();
	help();
}

void ChromosomeOneDigit::buildCodes() {
	// Can delete the codes
	canClean = true;

	// Make map
	codes = new map<char, char>();

	// Certain nucleotides
	codes->insert(map<char, char>::value_type('A', (char) 0));
	codes->insert(map<char, char>::value_type('C', (char) 1));
	codes->insert(map<char, char>::value_type('G', (char) 2));
	codes->insert(map<char, char>::value_type('T', (char) 3));

	// Common uncertain nucleotide
	// codes->insert(map<char, char>::value_type('N', (char) 4));

	// Uncertain nucleotides
	codes->insert(map<char, char>::value_type('R', codes->at('G')));
	codes->insert(map<char, char>::value_type('Y', codes->at('C')));
	codes->insert(map<char, char>::value_type('M', codes->at('A')));
	codes->insert(map<char, char>::value_type('K', codes->at('T')));
	codes->insert(map<char, char>::value_type('S', codes->at('G')));
	codes->insert(map<char, char>::value_type('W', codes->at('T')));
	codes->insert(map<char, char>::value_type('H', codes->at('C')));
	codes->insert(map<char, char>::value_type('B', codes->at('T')));
	codes->insert(map<char, char>::value_type('V', codes->at('A')));
	codes->insert(map<char, char>::value_type('D', codes->at('T')));
	codes->insert(map<char, char>::value_type('N', codes->at('C')));
	codes->insert(map<char, char>::value_type('X', codes->at('G')));
}

ChromosomeOneDigit::~ChromosomeOneDigit() {
	if (canClean) {
		codes->clear();
		delete codes;
	}
}

/**
 * This method converts nucleotides in the segments to single digit codes
 */
void ChromosomeOneDigit::encodeNucleotides() {

	for (int s = 0; s < segment->size(); s++) {
		int segStart = segment->at(s)->at(0);
		int segEnd = segment->at(s)->at(1);
		for (int i = segStart; i <= segEnd; i++) {

			if (codes->count(base[i]) > 0) {
				base[i] = codes->at(base[i]);
			} else {
				string msg = "Invalid nucleotide: ";
				msg.append(1, base[i]);
				throw InvalidInputException(msg);
			}
		}
	}

	// Digitize skipped segments
	int segNum = segment->size();
	if (segNum > 0) {
		// The first interval - before the first segment
		int segStart = 0;
		int segEnd = segment->at(0)->at(0) - 1;

		for (int s = 0; s <= segNum; s++) {
			for (int i = segStart; i <= segEnd; i++) {
				char c = base[i];

				if (c != 'N') {
					if (codes->count(c) > 0) {
						base[i] = codes->at(c);
					} else {
						string msg = "Invalid nucleotide: ";
						msg.append(1, c);
						throw InvalidInputException(msg);
					}
				}
			}

			// The regular intervals between two segments
			if (s < segNum - 1) {
				segStart = segment->at(s)->at(1) + 1;
				segEnd = segment->at(s + 1)->at(0) - 1;
			}
			// The last interval - after the last segment
			else if (s == segNum - 1) {
				segStart = segment->at(s)->at(1) + 1;
				segEnd = base.size() - 1;
			}
		}
	}
}

/**
 * Cannot be called on already finalized object.
 */
void ChromosomeOneDigit::makeR() {
	makeReverse();
	reverseSegments();
}

/**
 * Cannot be called on already finalized object.
 */
void ChromosomeOneDigit::makeRC() {
	makeComplement();
	makeReverse();
	reverseSegments();
}

void ChromosomeOneDigit::makeComplement() {
	map<char, char> complement;

	// Certain nucleotides
	complement.insert(map<char, char>::value_type((char) 0, (char) 3));
	complement.insert(map<char, char>::value_type((char) 1, (char) 2));
	complement.insert(map<char, char>::value_type((char) 2, (char) 1));
	complement.insert(map<char, char>::value_type((char) 3, (char) 0));

	// Unknown nucleotide
	complement.insert(map<char, char>::value_type('N', 'N'));

	// Convert a sequence to its complement
	int seqLen = base.size();
	for (int i = 0; i < seqLen; i++) {
		if (complement.count(base[i]) > 0) {
			base[i] = complement.at(base[i]);
		} else {
			cerr << "Error: The digit " << (char) base[i];
			cerr << " does not represent a base." << endl;
			exit(2);
		}
	}
}

void ChromosomeOneDigit::makeReverse() {
	int last = base.size() - 1;

	// Last index to be switched
	int middle = base.size() / 2;

	for (int i = 0; i < middle; i++) {
		char temp = base[last - i];
		base[last - i] = base[i];
		base[i] = temp;
	}
}

void ChromosomeOneDigit::reverseSegments() {
	int segNum = segment->size();
	int lastBase = size() - 1;

	// Calculate the coordinate on the main strand
	for (int i = 0; i < segNum; i++) {
		vector<int> * seg = segment->at(i);

		int s = lastBase - seg->at(1);
		int e = lastBase - seg->at(0);
		seg->clear();
		seg->push_back(s);
		seg->push_back(e);
	}

	// Reverse the regions within the list
	int lastRegion = segNum - 1;
	int middle = segNum / 2;
	for (int i = 0; i < middle; i++) {
		vector<int> * temp = segment->at(lastRegion - i);
		(*segment)[lastRegion - i] = segment->at(i);
		(*segment)[i] = temp;
	}
}

}
