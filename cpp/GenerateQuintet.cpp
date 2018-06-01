#include<iostream>
#include<algorithm>
#include<fstream>
 
#include "CSplitSys.h"
#include "auxilary.h"

using namespace std;

// Generate and order all 2|3 splits
void CSplitSys::GenerateQuintet()
{
     
     CQuintet aQuintet;
     
	 // push Binom(N,5) empty quintets to m_gQuintet
     for(unsigned long long n = 0; n != 10*Binom(m_iNumber,5); ++n)
             m_gQuintet.push_back(aQuintet);
             
     // List possible 5-tuples e.g. 54321.
     
     double dWei;
     
     unsigned int gLeft[2], gRight[3];
     unsigned int i1, i2, i3, i4, i5;
     for( i1 = m_iNumber; i1 != 4; --i1)
         for(i2 = i1-1; i2 != 3; --i2)
             for(i3 = i2-1; i3 != 2; --i3)
                 for(i4 = i3-1; i4 != 1; --i4)
                     for(i5 = i4-1; i5 != 0; --i5)
     {
           // Split of type 1: 21|543
           gLeft[0] = i4;
           gLeft[1] = i5;
           gRight[0] = i1;
           gRight[1] = i2;
           gRight[2] = i3;
           dWei = MinQuintet(gLeft, gRight);
            
           CQuintet tTemp1(gLeft, gRight, dWei);
           m_gQuintet[tTemp1.GetIndex()] = tTemp1;
           
           // Split of type 2: 31|542
           swap(gLeft[0], gRight[2]);
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp2(gLeft, gRight, dWei);
           m_gQuintet[tTemp2.GetIndex()] = tTemp2;
           
           // Split of type 3: 32|541
           swap(gLeft[1], gRight[2]);
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp3(gLeft, gRight, dWei);
           m_gQuintet[tTemp3.GetIndex()] = tTemp3;
           
           // Split of type 5: 42|531
           swap(gLeft[0], gRight[1]);
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp5(gLeft, gRight, dWei);
           m_gQuintet[tTemp5.GetIndex()] = tTemp5;
           
           // Split of type 4: 41|532
           swap(gLeft[1], gRight[2]);
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp4(gLeft, gRight, dWei);
           m_gQuintet[tTemp4.GetIndex()] = tTemp4;
           
           // Split of type 6: 43|521
           swap(gLeft[1], gRight[1]);
           swap(gRight[1], gRight[2]);
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp6(gLeft, gRight, dWei);
           m_gQuintet[tTemp6.GetIndex()] = tTemp6;
           
           // Split of type 9: 53|421
           swap(gLeft[0], gRight[0]);
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp9(gLeft, gRight, dWei);
           m_gQuintet[tTemp9.GetIndex()] = tTemp9;
           
           // Split of type 8: 52|431
           swap(gLeft[1], gRight[1]);
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp8(gLeft, gRight, dWei);
           m_gQuintet[tTemp8.GetIndex()] = tTemp8;
           
           // Split of type 7: 51|432
           swap(gLeft[1], gRight[2]);
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp7(gLeft, gRight, dWei);
           m_gQuintet[tTemp7.GetIndex()] = tTemp7;
           
          // Split of type 10: 54|321
           gLeft[1] = i2;
           gRight[0] = i3;
           gRight[1] = i4;
           gRight[2] = i5;
           dWei = MinQuintet(gLeft, gRight);   
           CQuintet tTemp10(gLeft, gRight, dWei);
           m_gQuintet[tTemp10.GetIndex()] = tTemp10;
     }
}
  
// EOF
