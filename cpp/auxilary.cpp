#include "auxilary.h"

typedef vector<string> vec;

// Return the number: n chooses k
unsigned long long Binom(unsigned int iN, unsigned int iK)
{
    if(iN < iK)
    {
          return 0;
    }
    else 
    {
         unsigned long long iBinom = 1;
         unsigned int i;
         
         for(i = 0; i < iK; ++i)
         {
             iBinom = iBinom * (iN - i);
         }
         
          for( i = 2; i < iK + 1; ++i)
          {
             iBinom = iBinom / i;
          } 
		  
          return iBinom;
    }    
}

//Remove spaces in the begin and end of a string  
string trim(string& s)  
{  
   const string &space =" \f\n\t\r\v" ;             // possible characters in the end of a string
   string r=s.erase(s.find_last_not_of(space)+1);   // remove the spaces in the end of the string
   return r.erase(0, r.find_first_not_of(space));   // remove the spaces in the begining of the string
}  

// Calculate the weights supporting s1|s2s3
unsigned int CalTriplet(const string& s1, const string& s2, const string& s3)
{
    unsigned int iNum=0;
    string::size_type i;
    
    for(i = 0; i != s1.size(); ++i)
    {
         if( ( s1[i] != s2[i] ) && ( s2[i] == s3[i]) ) iNum++;   //s1[i] \ne s2[i] = s3[i] indicates one support to s1|s2s3
    }
    
    return iNum;
}

// Calculate the weights supporting s1s2|s3s4
unsigned int CalQuartet(const string& s1, const string& s2, const string& s3, const string& s4)
{
    unsigned int iNum=0;
    string::size_type i;
    
    for( i = 0; i != s1.size(); ++i)
    {
         if(( s1[i] != s3[i] ) &&  ( s1[i] == s2[i] )  &&  ( s3[i] == s4[i] )) iNum++;  //s1[i] = s2[i] \ne s3[i]=s4[i] indicates one support to s1s2|s3s4
    }
    
    return iNum;
}
//Calculate the weights supporting s1s2s3|s4s5s6 
unsigned int CalSextet(const string& s1, const string& s2, const string& s3, const string& s4,const string& s5,const string& s6){
	
	unsigned int iNum=0;
    string::size_type i;
    
    for( i =0; i != s1.size(); ++i){
    	if(( s1[i] != s4[i] ) && ( s1[i]==s2[i] && s2[i]==s3[i]) && (s4[i]==s5[i] && s5[i]==s6[i])) iNum++;
	}
	
	return iNum;
}

bool cmp(unsigned int a,unsigned int b){
	return a>b;
}

//delete the same item in the vector vec
void vcomplete(vector<unsigned int>& vec,unsigned int temp1,unsigned int temp2) 
{
	vector<unsigned int>::iterator iter;
	for(iter = vec.begin(); iter != vec.end(); iter++){
		if( (temp1 == *iter) || (temp2 == *iter) ){
			vec.erase(iter);
		    --iter;
		} 
	}
}
//EOF
