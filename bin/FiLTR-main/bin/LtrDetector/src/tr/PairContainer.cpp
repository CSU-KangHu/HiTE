#include "PairContainer.h"
#include <iostream>
#include "Candidate.h"
#include <tuple>
#include <iterator>
#include <list>
#include <math.h>

using namespace std;

namespace tr{


    PairContainer::PairContainer(int minDistIn,int maxDistIn,int diffTolIn){

        minDistance = minDistIn;
        maxDistance = maxDistIn;
        diffTolerance = diffTolIn;

        allocate();

    }
    PairContainer::~PairContainer(){

        delete array;
    }

    void PairContainer::allocate(){

        containerSize = ((maxDistance-minDistance)/diffTolerance)+1;
       // cout<<"Size="<<containerSize<<endl;
       // cout<<"allocating!"<<endl;
        array = new std::list<Candidate*>[containerSize]; 
        
       // cout<<"allocate passed!"<<endl;
        
    }

    int PairContainer::computeIndex(int height){

        return ((height-minDistance)/diffTolerance);
    }
    /** Returns a candidate if a matching region (as defined by within height difference threshold)
 * is found. If no match is found, it hashes the candidate region and returns a null pointer.
 *
*/ 
    Candidate * PairContainer::hashOrReturn(Candidate * candidate){

         int height = candidate->getAbsHeight();
        // cout<<"height"<<height<<endl;
         int searchLoc = computeIndex(height);
         //cout<<"searchLoc="<<searchLoc<<endl;

        //cout<<"considering"<<endl;
        //candidate->printCandidate();
         
         std::list<Candidate *> midBin = array[searchLoc];

         if (candidate->getHeight() > 0)
         {
             midBin.emplace_front(candidate);
            // cout<<"emplacing"<<endl;
            // candidate->printCandidate();
             array[searchLoc] = midBin;
             return nullptr;
         }

         // cout<<"Emplacing"<<endl;
         

         else{

         if(!midBin.empty()){
                
                auto it = midBin.begin();

                while(it!=midBin.end())
               // for (auto it = midBin.begin(); it != midBin.end(); it++)
                {   //cout<<"Inside midBin"<<searchLoc<<endl;
                    Candidate *curr = *it;
                  //  cout<<"*it not deleted:"<<(curr!=nullptr)<<endl;

                    if ((abs(candidate->getHeight() + curr->getHeight()) <= diffTolerance) && ((curr->getStart() + curr->getHeight() + diffTolerance) >= candidate->getStart()) /*&& sameLength(candidate, curr)*/)
                    {   //cout<<1<<endl;
                        it = midBin.erase(it);
                        //cout<<2<<endl;
                        array[searchLoc] = midBin;
                        //cout << "returning" << endl;
                        //curr->printCandidate();
                        return curr;
                    }

                    

                    else{
                       // cout<<"3"<<endl;
                       // curr->printCandidate();
                        it = midBin.erase(it);
                        //cout <<"4"<<endl;
                        //cout<<"erasing"<<endl;
                       // array[searchLoc] = midBin;
                        //curr->printCandidate();
                    }

                    
                }
                array[searchLoc] = midBin;
              
        }

        if(searchLoc<containerSize -1){

            std::list<Candidate *> highBin = array[searchLoc+1];

            if(!highBin.empty()){

                //cout << "Inside highbin" << searchLoc+1 << endl;

                auto it = highBin.begin();

                while (it != highBin.end())
                //for (auto it = highBin.begin(); it != highBin.end(); it++)
                {
                    Candidate *curr = *it;

                    if ((abs(candidate->getHeight() + curr->getHeight()) <= diffTolerance) && ((curr->getStart() + curr->getHeight() + diffTolerance) >= candidate->getStart()) /*&& sameLength(candidate, curr)*/)
                    {
                       it = highBin.erase(it);
                        array[searchLoc+1] = highBin;
                        //cout << "returning" << endl;
                       // curr->printCandidate();


                        return curr;
                    }
                     
                    else {

                        it =highBin.erase(it);
                       // array[searchLoc+1] =highBin; 
                       // cout << "erasing" << endl;
                       // curr->printCandidate();
                    }
                    
                }
                array[searchLoc+1] = highBin;
            }
            
        }

        if(searchLoc >0){

            std::list<Candidate *> lowBin = array[searchLoc-1];

            if (!lowBin.empty())
            {   //cout<<"Inside lowBin"<<searchLoc-1<<endl;

                auto it = lowBin.begin();

                while (it != lowBin.end())
                    //for (auto it = lowBin.begin(); it != lowBin.end(); it++)
                    {
                        Candidate *curr = *it;

                        if ((abs(candidate->getHeight() + curr->getHeight()) <= diffTolerance) && ((curr->getStart() + curr->getHeight() + diffTolerance) >= candidate->getStart()) /*&& sameLength(candidate,curr)*/)
                        {
                            it = lowBin.erase(it);
                             array[searchLoc - 1] = lowBin;
                           // cout << "returning" << endl;
                           // curr->printCandidate();
                            return curr;
                        }

                        else{
                        it =lowBin.erase(it);
                        //array[searchLoc - 1] = lowBin;
                       // cout << "erasing" << endl;
                       // curr->printCandidate();
                    }
                    
                }

            array[searchLoc - 1] = lowBin;
            }
        }
        return nullptr;
        }

    
        
    }

    bool PairContainer :: sameLength(Candidate * curr, Candidate * next){

        int len1 = curr->getEnd()-curr->getStart();
        int len2 = next->getEnd()-next->getStart();

        int diff = abs(len2-len1);

        return diff<=20;


    }

    void PairContainer:: empty(){

        for(int i = 0;i<containerSize;i++){

            std::list<Candidate *> bin = array[i];
            if(!bin.empty()){
                for(auto it = bin.begin();it!=bin.end();it++){
                    Candidate * curr = *it;
                    cout<<"Start: "<< curr->getStart()<< "End: "<<curr->getEnd()<<"Height: "<<curr->getHeight()<<endl;
                }
            }
            cout<<"================================================================="<<endl;
        }
    }
}  
            
            
      
       