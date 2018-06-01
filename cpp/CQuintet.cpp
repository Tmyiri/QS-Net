#include <iostream> 
#include <algorithm>
#include "CQuintet.h"
#include "auxilary.h"

using namespace std;

// Constructor with parameters: m_gLblk = gLb; m_gRblk = gRb; m_dWeight= dWei
CQuintet::CQuintet(unsigned int gLb[2], unsigned int gRb[3], double dWei){
    m_gLblk[0] = gLb[0];
	m_gLblk[1] = gLb[1];
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];     
	m_gRblk[2] = gRb[2];                                                      
	m_dWeight = dWei;
	DecreaseOrder();
}

// copy constructor
CQuintet::CQuintet(const CQuintet& rhs){
	m_gLblk[0] = rhs.m_gLblk[0];
	m_gLblk[1] = rhs.m_gLblk[1];
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_gRblk[2] = rhs.m_gRblk[2];
	m_dWeight = rhs.m_dWeight;
}

// operator constructor
CQuintet& CQuintet::operator=(const CQuintet& rhs){
    m_gLblk[0] = rhs.m_gLblk[0];
	m_gLblk[1] = rhs.m_gLblk[1];
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_gRblk[2] = rhs.m_gRblk[2];
	m_dWeight = rhs.m_dWeight;
	return *this;
}

// set split label to be gLb|gRb 
void CQuintet::SetLabel(unsigned int gLb[2], unsigned int gRb[3]){
    m_gLblk[0] = gLb[0];
	m_gLblk[1] = gLb[1];
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];  
	m_gRblk[2] = gRb[2]; 
	DecreaseOrder();
}

// reorder quintet label according to the rules. e.g. 12|435 --> 21|543
void CQuintet::DecreaseOrder()
{  
   // reorder the left block 
   if (m_gLblk[0] < m_gLblk[1]) swap(m_gLblk[0], m_gLblk[1]);
   
   // reorder the right block
   if (m_gRblk[0] < m_gRblk[1]) swap(m_gRblk[0], m_gRblk[1]);
   
   if(m_gRblk[1] < m_gRblk[2]) swap(m_gRblk[1], m_gRblk[2]);
   
   if(m_gRblk[0] < m_gRblk[1]) swap(m_gRblk[0], m_gRblk[1]);
}

// calculate and return the index of quartet in the sereis 51|432 < 52|431 < 53|421 < 54|321 < 61|432 ...
unsigned int CQuintet::GetIndex() const
{  
	unsigned int index;
//	cout<<m_gLblk[0]<<' '<<m_gLblk[1]<<' '<<m_gRblk[0]<<' '<<m_gRblk[1]<<' '<<m_gRblk[2]<<endl;

    // Case: 21|543 
    if(m_gLblk[0] < m_gRblk[2]) index = 10*( Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gLblk[0]-1,2)+ Binom(m_gLblk[1]-1,1));
	
    // Case: 54|321    
	else if(m_gLblk[1] > m_gRblk[0]) index = 10*( Binom(m_gLblk[0]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[0]-1,3) + Binom(m_gRblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +9;
	
	// Case 3i|54j
    else if(m_gLblk[0] < m_gRblk[1])
   {
       // Case: 31|542 
       if(m_gLblk[1] < m_gRblk[2]) index = 10*( Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gLblk[0]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gLblk[1]-1,1)) +1;
	   
       // Case: 32|541
	   else  index = 10*( Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gLblk[0]-1,3) + Binom(m_gLblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +2;
   }
   
   // Case: 4i|5jk 
   else if(m_gLblk[0] < m_gRblk[0])
   {
       // Case: 41|532
       if(m_gLblk[1] < m_gRblk[2]) index = 10*( Binom(m_gRblk[0]-1,5) + Binom(m_gLblk[0]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gLblk[1]-1,1)) +3;
	   
       // Case: 42|531
	   else if(m_gLblk[1] < m_gRblk[1])   index = 10*( Binom(m_gRblk[0]-1,5) + Binom(m_gLblk[0]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gLblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +4;
	   
       // Case: 43|521
	   else   index = 10*( Binom(m_gRblk[0]-1,5) + Binom(m_gLblk[0]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gRblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +5;
   }
   
   // 5i|4jk
   else
   {
       // Case: 51|432
       if(m_gLblk[1] < m_gRblk[2]) index = 10*( Binom(m_gLblk[0]-1,5) + Binom(m_gRblk[0]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gLblk[1]-1,1)) +6;
	   
       // Case: 52|431
	   else if(m_gLblk[1] < m_gRblk[1])  index = 10*( Binom(m_gLblk[0]-1,5) + Binom(m_gRblk[0]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gLblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +7;
		  
	    // Case: 53|421   
	   else index = 10*( Binom(m_gLblk[0]-1,5) + Binom(m_gRblk[0]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gRblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +8;
   }
   
   return index;
}   

//EOF

