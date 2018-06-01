#include <algorithm>
#include "CQuartet.h"
#include "auxilary.h"

using namespace std;

// Constructor with parameters: m_gLblk = gLb; m_gRblk = gRb; m_dWeight= dWei
CQuartet::CQuartet(unsigned int gLb[2], unsigned int gRb[2], double dWei){
    m_gLblk[0] = gLb[0];
	m_gLblk[1] = gLb[1];
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];                                                          
	m_dWeight = dWei;
	DecreaseOrder();
}

// Copy constructor
CQuartet::CQuartet(const CQuartet& rhs){
	m_gLblk[0] = rhs.m_gLblk[0];
	m_gLblk[1] = rhs.m_gLblk[1];
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_dWeight = rhs.m_dWeight;
}

// Operator constructor
CQuartet& CQuartet::operator=(const CQuartet& rhs){
    m_gLblk[0] = rhs.m_gLblk[0];
	m_gLblk[1] = rhs.m_gLblk[1];
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_dWeight = rhs.m_dWeight;
	return *this;
}
 
// Set split label to be gLb|gRb 
void CQuartet::SetLabel(unsigned int gLb[2], unsigned int gRb[2]){
    m_gLblk[0] = gLb[0];
	m_gLblk[1] = gLb[1];
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];  
	DecreaseOrder();
}

// Reorder quartet label according to the format 12|34 = 43|21
void CQuartet::DecreaseOrder()
{   
   if (m_gLblk[0] < m_gLblk[1]) swap(m_gLblk[0], m_gLblk[1]);
   
   if (m_gRblk[0] < m_gRblk[1]) swap(m_gRblk[0], m_gRblk[1]);
   
   if(m_gLblk[0]<m_gRblk[0])
   {
        swap(m_gLblk[0], m_gRblk[0]);
        swap(m_gLblk[1], m_gRblk[1]);
   }
}

// Calculate and return the index of quartet in the sereis 41|32 < 42|31 < 43|21 < 51|32 ...
unsigned int CQuartet::GetIndex() const
{   
	unsigned int index;
	
    // Case: 43|21
    if(m_gLblk[1] > m_gRblk[0])
        index = 3* (Binom(m_gLblk[0]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gRblk[0]-1,2) + Binom(m_gRblk[1]-1,1))+2;
		
	// Case: 42|31
    else if (m_gLblk[1] > m_gRblk[1])
		index =  3* (Binom(m_gLblk[0]-1,4) + Binom(m_gRblk[0]-1,3) + Binom(m_gLblk[1]-1,2) + Binom(m_gRblk[1]-1,1))+1;
	
	 // Case: 41|32
	else
		index =  3* (Binom(m_gLblk[0]-1,4) + Binom(m_gRblk[0]-1,3) + Binom(m_gRblk[1]-1,2) + Binom(m_gLblk[1]-1,1));
	
	return index;
}   

//EOF

