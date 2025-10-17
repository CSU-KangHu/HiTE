#ifndef PAIRCONTAINER_H_
#define PAIRCONTAINER_H_
#include "Candidate.h"
#include <vector>
#include<list>

using namespace std;

namespace tr{

    class PairContainer{

        private:
            std::list<Candidate*> * array;
            int maxDistance;
            int minDistance;
            int diffTolerance;
            bool sameLength( Candidate *, Candidate *);
            int containerSize;
            void allocate();
            int computeIndex(int);

        public:
            PairContainer(int,int,int);
            virtual ~PairContainer();
            Candidate * hashOrReturn(Candidate *);
            void empty();
           // void unHash(int);
           // Container removeMatch();
           // tuple<int,int> at(int);


    };

}
#endif