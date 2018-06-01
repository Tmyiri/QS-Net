#include <iostream> 
#include <algorithm>
#include "CSeptet.h"
#include "auxilary.h"

using namespace std;

// Constructor with parameters: m_gLblk = gLb; m_gRblk = gRb; m_dWeight= dWei
CSeptet::CSeptet(unsigned int gLb[3], unsigned int gRb[4], double dWei){
    m_gLblk[0] = gLb[0];
	m_gLblk[1] = gLb[1];
	m_gLblk[2] = gLb[2];
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];     
	m_gRblk[2] = gRb[2]; 
	m_gRblk[3] = gRb[3];                                                     
	m_dWeight = dWei;
	DecreaseOrder();
}

// copy constructor
CSeptet::CSeptet(const CSeptet& rhs){
	m_gLblk[0] = rhs.m_gLblk[0];
	m_gLblk[1] = rhs.m_gLblk[1];
	m_gLblk[2] = rhs.m_gLblk[2];
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_gRblk[2] = rhs.m_gRblk[2];
	m_gRblk[3] = rhs.m_gRblk[3];
	m_dWeight = rhs.m_dWeight;
}

// operator constructor
CSeptet& CSeptet::operator=(const CSeptet& rhs){
    m_gLblk[0] = rhs.m_gLblk[0];
	m_gLblk[1] = rhs.m_gLblk[1];
	m_gLblk[2] = rhs.m_gLblk[2];
	m_gRblk[0] = rhs.m_gRblk[0];
	m_gRblk[1] = rhs.m_gRblk[1];
	m_gRblk[2] = rhs.m_gRblk[2];
	m_gRblk[3] = rhs.m_gRblk[3];
	m_dWeight = rhs.m_dWeight;
	return *this;
}

// set split label to be gLb|gRb 
void CSeptet::SetLabel(unsigned int gLb[3], unsigned int gRb[4]){
    m_gLblk[0] = gLb[0];
	m_gLblk[1] = gLb[1];
	m_gLblk[2] = gLb[2];
    m_gRblk[0] = gRb[0];
	m_gRblk[1] = gRb[1];  
	m_gRblk[2] = gRb[2]; 
	m_gRblk[3] = gRb[3]; 
	DecreaseOrder();
}

// reorder septet label according to the rules. e.g. 123|4567 --> 321|7654
void CSeptet::DecreaseOrder()
{  
   sort(m_gLblk,m_gLblk+3,cmp);
   
   sort(m_gRblk,m_gRblk+4,cmp);
}

// calculate and return the index of septet in the sereis 321|7654 < 421|7653 < 431|7652 < 432|7651 < 521|7643...
unsigned int CSeptet::GetIndex() const
{  
    unsigned int index;
	
    // Case: 321|7654 
    if(m_gLblk[0] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gRblk[2]-1,5) + Binom(m_gRblk[3]-1,4) + Binom(m_gLblk[0]-1,3) + Binom(m_gLblk[1]-1,2)+ Binom(m_gLblk[2]-1,1));
	
    // Case: 765|4321    
	else if(m_gLblk[2] > m_gRblk[0]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gLblk[1]-1,6) +Binom(m_gLblk[2]-1,5) + Binom(m_gRblk[0]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1)) +34;
	
	// Case 4ij|765k
    else if(m_gLblk[0] < m_gRblk[2])
   {
   	   // Case: 421|7653
       if(m_gLblk[1] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gRblk[2]-1,5) + Binom(m_gLblk[0]-1,4) + Binom(m_gRblk[3]-1,3) + Binom(m_gLblk[1]-1,2)+ Binom(m_gLblk[2]-1,1))+1;
       
       //Case: 431|7652
       else if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gRblk[2]-1,5) + Binom(m_gLblk[0]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+2;
       
       //Case: 432|7651
	   else    index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gRblk[2]-1,5) + Binom(m_gLblk[0]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+3; 
   }
   
   // Case: 5ij|76jk
   else if(m_gLblk[0] < m_gRblk[1])
   {
       //Case: 521|7643
       if(m_gLblk[1] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gLblk[0]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gRblk[3]-1,3) + Binom(m_gLblk[1]-1,2)+ Binom(m_gLblk[2]-1,1))+4;
       
       //Case: 53i|764j
       else if(m_gLblk[1] < m_gRblk[2]){
       	    //Case: 531|7642
       	    if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gLblk[0]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+5;
       	    
       	    //Case: 532|7641
       	    else    index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gLblk[0]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+6;    	   
	   }
	   
	   //Case: 54i|76jk
	   else{
	   	   //Case:541|7632
	   	   if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gLblk[0]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+7;
	   	   
	   	   //Case: 542|7631
	   	   else if(m_gLblk[2] < m_gRblk[2]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gLblk[0]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+8;
	   
	       //Case: 543|7621
	       else     index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gRblk[1]-1,6) +Binom(m_gLblk[0]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gLblk[2]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+9;
	   }
   }
   
   // Case: 6ij|7jkl
   else if(m_gLblk[0] < m_gRblk[0])
   {
       //Case: 621|7543
       if(m_gLblk[1] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gRblk[3]-1,3) + Binom(m_gLblk[1]-1,2)+ Binom(m_gLblk[2]-1,1))+10;
       
       //Case: 63i|754j
       else if(m_gLblk[1] < m_gRblk[2])
	   {
       	    //Case: 631|7542
       	    if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+11;
       	    
       	    //Case: 632|7541
       	    else    index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+12;    	   
	   }
	   
	   //Case: 64i|75jk
	   else if(m_gLblk[1] < m_gRblk[1])
	   {
	   	   //Case:641|7532
	   	   if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+13;
	   	   
	   	   //Case: 642|7531
	   	   else if(m_gLblk[2] < m_gRblk[2]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+14;
	   
	       //Case: 643|7521
	       else index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gLblk[2]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+15;
	   	return index;
	   }
	   
	   //Case: 65i|7jkl
	   else
	   {
	   	   //Case:651|7432
	   	   if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+16;
	   	   
	   	   //Case: 652|7431
	   	   else if(m_gLblk[2] < m_gRblk[2]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+17;
	   	   
	   	   //Case: 653|7421
	   	   else if(m_gLblk[2] < m_gRblk[1]) index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gLblk[2]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+18;
	   	   
	   	   //Case: 654|7321
	   	   else     index = 35*( Binom(m_gRblk[0]-1,7) +Binom(m_gLblk[0]-1,6) +Binom(m_gLblk[1]-1,5) + Binom(m_gLblk[2]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+19;
	   }
   }
   
   //Case: 7ij|klmn
   else
   {
   	   //Case: 721|6543
       if(m_gLblk[1] < m_gRblk[3]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gRblk[3]-1,3) + Binom(m_gLblk[1]-1,2)+ Binom(m_gLblk[2]-1,1))+20;
       
       //Case: 73i|754j
       else if(m_gLblk[1] < m_gRblk[2])
	   {
       	    //Case: 731|6542
       	    if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+21;
       	    
       	    //Case: 732|6541
       	    else    index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gRblk[2]-1,4) + Binom(m_gLblk[1]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+22;  	   
	   }
	   
	   //Case: 74i|65jk
	   else if(m_gLblk[1] < m_gRblk[1])
	   {
	   	   //Case:741|6532
	   	   if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+23;
	   	   
	   	   //Case: 742|6531
	   	   else if(m_gLblk[2] < m_gRblk[2]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+24;
	   
	       //Case: 743|6521
	       else     index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gRblk[1]-1,5) + Binom(m_gLblk[1]-1,4) + Binom(m_gLblk[2]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+25;
	   }
	   
	   //Case: 75i|6jkl
	   else if(m_gLblk[1] < m_gRblk[0])
	   {
	   	   //Case: 751|6432
	   	   if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+26;
	   	   
	   	   //Case: 752|6431
	   	   else if(m_gLblk[2] < m_gRblk[2]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+27;
	   	   
	   	   //Case: 753|6421
	   	   else if(m_gLblk[2] < m_gRblk[1]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gLblk[1]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gLblk[2]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+28;
	   	   
	   	   //Case: 754|6321
	   	   else     index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gRblk[0]-1,6) +Binom(m_gLblk[1]-1,5) + Binom(m_gLblk[2]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+29;
	   }
	   
	   //Case: 76i|5jkl
	   else
	   {
	   	   //Case: 761|5432
	   	   if(m_gLblk[2] < m_gRblk[3]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gLblk[1]-1,6) +Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gRblk[3]-1,2)+ Binom(m_gLblk[2]-1,1))+30;
	   	   
	   	   //Case: 762|5431
	   	   else if(m_gLblk[2] < m_gRblk[2]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gLblk[1]-1,6) +Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gRblk[2]-1,3) + Binom(m_gLblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+31;
	   	   
	   	   //Case: 763|5421
	   	   else if(m_gLblk[2] < m_gRblk[1]) index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gLblk[1]-1,6) +Binom(m_gRblk[0]-1,5) + Binom(m_gRblk[1]-1,4) + Binom(m_gLblk[2]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+32;
	   	   
	   	   //Case: 764|5321
	   	   else     index = 35*( Binom(m_gLblk[0]-1,7) +Binom(m_gLblk[1]-1,6) +Binom(m_gRblk[0]-1,5) + Binom(m_gLblk[2]-1,4) + Binom(m_gRblk[1]-1,3) + Binom(m_gRblk[2]-1,2)+ Binom(m_gRblk[3]-1,1))+33;
	   }
   }
   
   return index;
}   

//EOF
