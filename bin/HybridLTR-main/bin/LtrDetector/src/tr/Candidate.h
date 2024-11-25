#ifndef CANDIDATE_H_
#define CANDIDATE_H_
using namespace std;

namespace tr{

class Candidate{

    private:
        int start;
        int end;
        int height;
    public:
        int getStart();
        int getEnd();
        int getHeight();
        int getAbsHeight();
        Candidate(int,int,int);
        ~Candidate();
        void printCandidate();

};
}
#endif