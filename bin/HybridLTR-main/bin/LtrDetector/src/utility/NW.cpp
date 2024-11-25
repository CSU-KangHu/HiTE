/* -*- C++ -*-
 *
 * needleman_wunsch.cpp
 *
 * Author: Benjamin T James
 */
#include "NW.h"
#include <iostream>

//flags that can be combined
#define HORIZ 1
#define VERT  2
#define DIAG  4
using namespace std;
void needleman_wunsch::fill(int i, int j)
{
	if (i == 0 || j == 0) {
		if (i == j) {
			int offset = at(i, j);
			score[offset] = 0;
			direction[offset] = DIAG; // for backtracking
			horiz_gap_len[offset] = 0;
			vert_gap_len[offset] = 0;
		} else if (i == 0) {
			int offset = at(0, j);
			int last_offset = at(0, j-1);
			score[offset] = score[last_offset] + gap(j);
			horiz_gap_len[offset] = 0;
			vert_gap_len[offset] = j;
			direction[offset] = VERT;
		} else { // j == 0
			int offset = at(i, 0);
			int last_offset = at(i-1, 0);
			score[offset] = score[last_offset] + gap(i);
			horiz_gap_len[offset] = i;
			vert_gap_len[offset] = 0;
			direction[offset] = HORIZ;
		}
		return;
	}
	int i_diag = at(i-1, j-1);
	int i_horiz = at(i-1, j);
	int i_vert = at(i, j-1);
	int i_cur = at(i, j);

	int hlen = horiz_gap_len[i_horiz] + 1;
	int vlen = vert_gap_len[i_vert] + 1;

	int diag_score = score[i_diag] + match_score(s1[i], s2[j]);
	int horiz_score = score[i_horiz] + gap(hlen);
	int vert_score = score[i_vert] + gap(vlen);
	score[i_cur] = std::max(std::max(diag_score, horiz_score), vert_score);
	direction[i_cur] = 0;

	// we could match multiple high scores
	if (score[i_cur] == diag_score) {
		direction[i_cur] |= DIAG;
	}
	if (score[i_cur] == vert_score) {
		direction[i_cur] |= VERT;
		vert_gap_len[i_cur] = vlen;
	} else {
		vert_gap_len[i_cur] = 0;
	}
	if (score[i_cur] == horiz_score) {
		direction[i_cur] |= HORIZ;
		horiz_gap_len[i_cur] = hlen;
	} else {
		horiz_gap_len[i_cur] = 0;
	}
}

std::pair<std::string, std::string>
needleman_wunsch::backtrack()
{
	std::string a1 = "", a2 = "";
	int cur_i = l1 - 1;
	int cur_j = l2 - 1;
	while (cur_i >= 0 && cur_j >= 0) {
		uint8_t dir = direction[at(cur_i, cur_j)];
		if (dir & DIAG) {
			a1 += s1[cur_i--];
			a2 += s2[cur_j--];
		} else if (dir & HORIZ) {
			a1 += s1[cur_i--];
			a2 += '-';
		} else if (dir & VERT) {
			a1 += '-';
			a2 += s2[cur_j--];
		}
	}
	std::string r1(a1.rbegin(), a1.rend());
	std::string r2(a2.rbegin(), a2.rend());
	return std::make_pair(r1, r2);
}


std::pair<std::string, std::string>
needleman_wunsch::align()
{
	for (int i = 0; i < l1; i++) {
		for (int j = 0; j < l2; j++) {
			fill(i, j);
		}
	}
	return backtrack();
}
double needleman_wunsch::identity(std::pair<std::string, std::string> alignment) const
{
	int len = alignment.first.length();
	double count = 0;
	for (int i = 0; i < len; i++) {
		if (alignment.first[i] == alignment.second[i]) {
			count++;
		}
	}
	return (double)count / len;
}

int needleman_wunsch::gap(int gaplen) const
{
	return sigma + (gaplen - 1) * epsilon;
}

int transform1(char a)
{
	int ret = -1;
	switch (a) {
	case 'a':
	case 'A':
		ret = 0;
		break;
	case 'c':
	case 'C':
		ret = 1;
		break;
	case 'g':
	case 'G':
		ret = 2;
		break;
	case 't':
	case 'T':
		ret = 3;
		break;
	default:
		ret = a;
		break;
	}
	return ret;
}
int needleman_wunsch::match_score(char a, char b) const
{
        int ta = transform1(a);
	int tb = transform1(b);
	return scoring_matrix[ta][tb];
}
needleman_wunsch::needleman_wunsch(const std::string &s1_, const std::string& s2_, const int m[4][4], int sigma_, int epsilon_)
{
	int l1_ = s1_.length();
	int l2_ = s2_.length();
	if (l1_ >= l2_) {
		l1 = l1_;
		l2 = l2_;
		s1 = s1_;
		s2 = s2_;
	} else {
		l1 = l2_;
		l2 = l1_;
		s1 = s2_;
		s2 = s1_;
	}
	sigma = sigma_;
	epsilon = epsilon_;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			scoring_matrix[i][j] = (int)m[i][j];
		}
	}
	int matlen = l1 * l2;
	score = new int[matlen];
	direction = new uint8_t[matlen];
	horiz_gap_len = new int[matlen];
	vert_gap_len = new int[matlen];
}

int main(){
   // const std::string string1 = "GCTCCCTGGAGGTAGGAGCGTACGTCGGAGGCGTCAGCAGCATGTCTCAGTGGGACCGTTGGTTCCATTGTTAGATTGTTGATCACTACTTGTGTGTTTTGAAGCTATCTAACGCACGTTGACTTAGTATCTCTATAGTTCCGATTGACTAATTACACCTCGAGTACATTTAAGTGACTCTTAGGTAATGCGTTAGGCAAGCAAAATCTGACGCCCACGTACACCGATGCCCATAGAGTCAGGAGGGGCATTAGTCTCGGACGACATCGACGGCGATATCAGGCTTTCTACTCCGCCCTTAAGGACACCGAAACTCGGTTATGAAGAACGCGTGGCATTAAGCCGCTGCCCACTGTGGTTGTCAGGCTCGTGTATCAAGCATAACATAACAGCGGGACCAAGCATTCAGGCGTTTTGATTAAGACCGATGTACCAAGAACGACGAGGTACGGGGTGACAACAAAGTTCTCTAAGGATACATGATTGGGGGCTCAGCAATGAATCTGATCTTCCATAGAAGGATAGTACCTCTCCGTAGTCTCACTTCGCGGACTGCCGTTCAGTTTTCCTATACATTGCTCTCGAATTGCGCGTTTAAGTTTGCTTCAGTTGGGAACACGATTTTGGTGTAGAACGTTAGAAAAGTAACTCAGAGGGGTGCGGTGTAAGTTGTTCACCTTCTGCTGGGCAATCACGGTGAGCCCTTCCAGCGTGCCACGAATTCGATACCCCACGTGATCTAGCTGGCTGGCCCAACCGCATGTTGGAACGTGAGACGGCCAGACACCGAGCACAGGTATTGACCTCCGGGCAAACACTCGGATCGATCTTCGTACAACGTCTTTGTGTTTCCCTATTGAATTTTCCCCGCGTCATGTTCGATCCATCACGACCAACGAGGTGGACCAAGGAGTGAATTCTGAAGATCCGAAACTTTTTAATGTAAACTACCGATGTGAAAAACCAAAAAATTCGTTAGGCTTACTACCAGAATAGAGTT";
   // const std::string string2 = "cCtCtCCtGGAGGTAGGAGCGTACGaCGGAGGCGTCAGCAGCATacGTCgCAGTGGGACccGTTGGgTCCATTGTTAGATTGTTGATCACTACTTGaTgtTTTGagCTATCTAACGccGTTGACTTAGTtATCTCTATAGtttcGATTGAaTAATTaaCAtCTCgAGttActTTAAGTGACTCTTAGGTAATGCGTtTAGGCaAGCAAAATCTGACGCCCACGTtccGATGCCCATAGgAGTCAGcGAggGCATTAGTCTCGGACGaaTCGcACGGCGATATCAGGCTTtttACTCgCGCCCTTAAGGACACCGAAACTCgGTttATGAAGtACgcgcGTcggATTAacaGCTGCccCACttGTGGTTGTCAGaGCTCGTGTATCAAggCATAACATaACAGCGGGACCAAGCtgcATTCAGaGCGaTTtgATTAAGACcaTGTACCAaaGgAACGACGcAGgacGGGGTGaaACAAatCTCTAAGGATACATGATtGGGGGcTCcAGCAATGAATCTgATCTTCCATAGAAaGGATAGtcCTCTCCGcAGTCTCatTCGCgaCTGCCGTTcatTTTCCTATcgACATtcTCTcaATTGaCGCGTTTagTTTGcCTTCAGTgGGGAACACGATTTTGGTGTAGAACGTTAGAAAAGTAAccAGAGGGGTGCGGTGTgAgtGTTCACCTTCTGCTGGGCAATCgCGGTaGAGCCCTTCCAGCGTGCCACGAatgATACccCCAcggTGATCTAGCTGGCTGGCCCAACCGCAgGTTGGAActGAGACGGCCAGACACccGAGCacacAGGTATTGACCTCCGGGCAAACACTcgATCacggATCTTCGTACAACGTCTTTtgGTTTCCCTATTggAATTTTccCGCGTCATGTTCGATCCATCACGACCAACGAGGaTGgACCAAGGAGcgaTTCTGAAGATCcgaCgaaaACTTttTTAATGTAAACTACCGATGTGAAAAACCAAAatTCGTTAGGCTcTACTACCAGAATAGAGTT";
    string string1 = "GCTCCCTGGAGGTAGG";
    string string2 ="cCtCtCCtGGAGGTAG"; 

    const std::string& in1 = string1;
    const std::string& in2 = string2;
     
      const int array [4][4]={{1,-1,-1,1},{-1,1,-1,-1},{-1,-1,1,-1},{-1,-1,-1,1}};
 

     needleman_wunsch * alignment = new needleman_wunsch(in1,in2,array,4,1);
   cout << "ALIGNMENT IS" << alignment->align().first<<endl;
     std::cout  << alignment->align().second<<std::endl;
}