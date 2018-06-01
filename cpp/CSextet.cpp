#include <algorithm>
#include "CSextet.h"
#include "auxilary.h"

using namespace std;

// Constructor with parameters: m_gLblk = gLb; m_gRblk = gRb; m_dWeight= dWei
CSextet::CSextet(unsigned int gLb[3], unsigned int gRb[3], double dWei){
    m_gLblk[0] = gLb[0];
	m_gLblk[1] = gLb[1];
	m_gLblk[2] = gLb[2];
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];  
	m_gRblk[2] = gRb[2];                                                        
	m_dWeight = dWei;
	DecreaseOrder();
}

// Copy constructor
CSextet::CSextet(const CSextet& rhs){
	m_gLblk[0] = rhs.m_gLblk[0];
	m_gLblk[1] = rhs.m_gLblk[1];
	m_gLblk[2] = rhs.m_gLblk[2];
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_gRblk[2] = rhs.m_gRblk[2];
	m_dWeight = rhs.m_dWeight;
}

// Operator constructor
CSextet& CSextet::operator=(const CSextet& rhs){
    m_gLblk[0] = rhs.m_gLblk[0];
	m_gLblk[1] = rhs.m_gLblk[1];
	m_gLblk[2] = rhs.m_gLblk[2];
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_gRblk[2] = rhs.m_gRblk[2];
	m_dWeight = rhs.m_dWeight;
	return *this;
}
 
// Set split label to be gLb|gRb 
void CSextet::SetLabel(unsigned int gLb[3], unsigned int gRb[3]){
    m_gLblk[0] = gLb[0];
	m_gLblk[1] = gLb[1];
	m_gLblk[2] = gLb[2];
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1]; 
	m_gRblk[2] = gRb[2]; 
	DecreaseOrder();
}

// Reorder sextet label according to the format 123|456 --> 654|321
void CSextet::DecreaseOrder()
{   
   if (m_gLblk[0] < m_gLblk[1]) swap(m_gLblk[0], m_gLblk[1]);
   if (m_gLblk[0] < m_gLblk[2]) swap(m_gLblk[0], m_gLblk[2]);
   if (m_gLblk[1] < m_gLblk[2]) swap(m_gLblk[1], m_gLblk[2]);
        
   if (m_gRblk[0] < m_gRblk[1]) swap(m_gRblk[0], m_gRblk[1]);
   if (m_gRblk[0] < m_gRblk[2]) swap(m_gRblk[0], m_gRblk[2]);
   if (m_gRblk[1] < m_gRblk[2]) swap(m_gRblk[1], m_gRblk[2]);
   
   if(m_gLblk[0]<m_gRblk[0])
   {
        swap(m_gLblk[0], m_gRblk[0]);
        swap(m_gLblk[1], m_gRblk[1]);
        swap(m_gLblk[2], m_gRblk[2]);
   }
}

// calculate and return the index of sextet in the sereis 621|543 < 631|542 < 632|541 < 641|532 ...
unsigned int CSextet::GetIndex() const
{   
	unsigned int index;
	
   // Case: 621|543 
    if(m_gLblk[1] < m_gRblk[2]) index = 10*( Binom(m_gLblk[0]-1,6) + Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gRblk[2]-1,3)+ Binom(m_gLblk[1]-1,2)+ Binom(m_gLblk[2]-1,1));
	
    // Case: 654|321    
	else if(m_gLblk[2] > m_gRblk[0]) index = 10*( Binom(m_gLblk[0]-1,6) + Binom(m_gLblk[1]-1,5) + Binom(m_gLblk[2]-1,4) + Binom(m_gRblk[0]-1,3)+ Binom(m_gRblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +9;
	
	// Case£º 63i|54j
    else if(m_gLblk[1] < m_gRblk[1])
   {
       // Case: 631|542 
       if(m_gLblk[2] < m_gRblk[2]) index = 10*( Binom(m_gLblk[0]-1,6) + Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gLblk[1]-1,3)+ Binom(m_gRblk[2]-1,2)+  Binom(m_gLblk[2]-1,1)) +1;
	   
       // Case: 632|541
	   else  index = 10*( Binom(m_gLblk[0]-1,6)+ Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[2]-1,1)) +2;
   }
   
   // Case: 64i|5jk 
   else if(m_gLblk[1] < m_gRblk[0])
   {
       // Case: 641|532
       if(m_gLblk[2] < m_gRblk[2]) index = 10*( Binom(m_gLblk[0]-1,6)+ Binom(m_gRblk[0]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gLblk[2]-1,1)) +3;
	   
       // Case: 642|531
	   else if(m_gLblk[2] < m_gRblk[1])   index = 10*( Binom(m_gLblk[0]-1,6)+ Binom(m_gRblk[0]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[2]-1,1)) +4;
	   
       // Case: 643|521
	   else   index = 10*( Binom(m_gLblk[0]-1,6)+ Binom(m_gRblk[0]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gLblk[2]-1,3) + Binom(m_gRblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +5;
   }
   
   // Case: 65i|4jk
   else
   {
       // Case: 651|432
       if(m_gLblk[2] < m_gRblk[2]) index = 10*( Binom(m_gLblk[0]-1,6)+ Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[0]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gLblk[2]-1,1)) +6;
	   
       // Case: 652|431
	   else if(m_gLblk[2] < m_gRblk[1])  index = 10*( Binom(m_gLblk[0]-1,6)+ Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[0]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[2]-1,1)) +7;
		  
	    // Case: 653|421   
	   else index = 10*( Binom(m_gLblk[0]-1,6)+ Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[0]-1,4) + Binom(m_gLblk[2]-1,3) + Binom(m_gRblk[1]-1,2)+ Binom(m_gRblk[2]-1,1)) +8;
   }
   	
	return index;
}   

//EOF

