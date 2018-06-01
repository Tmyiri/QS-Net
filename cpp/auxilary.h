#ifndef _AUXILARY_H
#define _AUXILARY_H

//Body of auxilary.h
#include <vector>
#include <string>

using namespace std;
                    
extern unsigned long long Binom(unsigned int, unsigned int);                                 // Return the number: n chooses m
extern string trim(string&);                                                                 // Remove spaces in the begining and end of a string
extern unsigned int CalTriplet(const string&, const string&, const string&);                 // Calculate the weights supporting s1|s2s3
extern unsigned int CalQuartet(const string&, const string&, const string&, const string&);  // Calculate the weights supporting s1s2|s3s4
extern unsigned int CalSextet(const string&, const string&, const string&, const string&, const string&, const string&);  // Calculate the weights supporting s1s2s3|s4s5s6
extern bool cmp(unsigned int, unsigned int);
extern void vcomplete(vector<unsigned int>&, unsigned int, unsigned int);                    // delete the same item in the vector vec

// Check if all the elements in T are different
template <typename T> bool DiffVec(const vector<T>& gVec)
{
    typename vector<T>::const_iterator it1, it2;
    
    for(it1 = gVec.begin(); it1 != gVec.end()-1; ++it1)
        for(it2 = it1+1; it2 != gVec.end(); ++it2)
       {
            if(*it1 == *it2) return false;
        }
        
   return true;   
}

// Check if all the elements in T are the same
template <typename T> bool SameVec(const vector<T>& gVec)
{
    typename vector<T>::const_iterator it1, it2;
    
    for(it1 = gVec.begin(); it1 != gVec.end()-1; ++it1)
        for(it2 = it1+1; it2 != gVec.end(); ++it2)
       {
            if(*it1 != *it2) return false;
        }
        
   return true;   
}

#endif
//EOF
