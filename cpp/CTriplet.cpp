#include <algorithm>
#include "CTriplet.h"
#include "auxilary.h"

using namespace std;

// Constructor with parameters: m_gLblk = gLb; m_gRblk = gRb; m_dWeight= dWei
CTriplet::CTriplet(unsigned int gLb, unsigned int gRb[2], double dWei){
    m_gLblk = gLb;
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];                                                          
	m_dWeight = dWei;
	DecreaseOrder();
}

// Copy constructor
CTriplet::CTriplet(const CTriplet& rhs){
	m_gLblk = rhs.m_gLblk;
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_dWeight = rhs.m_dWeight;
}

// Operator constructor
CTriplet& CTriplet::operator=(const CTriplet& rhs){
    m_gLblk = rhs.m_gLblk;
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_dWeight = rhs.m_dWeight;
	return *this;
}
 
// Set split label to be gLb|gRb 
void CTriplet::SetLabel(unsigned int gLb, unsigned int gRb[2]){
    m_gLblk = gLb;
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];  
	DecreaseOrder();
}


// Reorder quartet label according to the format 1|34 = 1|43
void CTriplet::DecreaseOrder()
{   
   if (m_gRblk[0] < m_gRblk[1]) swap(m_gRblk[0], m_gRblk[1]);
}

// Calculate and return the index of quartet in the sereis 41|32 < 42|31 < 43|21 < 51|32 ...
unsigned int CTriplet::GetIndex() const
{   
	unsigned int index; 
	
    // Case: 3|21
    if(m_gLblk > m_gRblk[0]) index = 3* (Binom(m_gLblk-1,3) + Binom(m_gRblk[0]-1,2) + Binom(m_gRblk[1]-1,1))+2;
	
    // Case: 2|31
    else if (m_gLblk> m_gRblk[1]) index = 3* (Binom(m_gRblk[0]-1,3) + Binom(m_gLblk-1,2) + Binom(m_gRblk[1]-1,1))+1;
	
    // Case: 1|32
	else  index = 3* ( Binom(m_gRblk[0]-1,3) + Binom(m_gRblk[1]-1,2) + Binom(m_gLblk-1,1));

    return index;
}   

//EOF

