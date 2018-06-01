#include<iostream>
#include<algorithm>
#include<fstream>
 
#include "CSplitSys.h"
#include "auxilary.h"

using namespace std;


// Generate and order all 3|4 splits
void CSplitSys::GenerateSeptet()
{
     
     CSeptet sSeptet;
     
	 // push Binom(N,7) empty septets to m_gSeptet
     for(unsigned long long n = 0; n != 35*Binom(m_iNumber,7); ++n)
             m_gSeptet.push_back(sSeptet);
             
     // List possible 7-tuples e.g. 1234567.
     
     double dWei;
     
     unsigned int SepgLeft[3], SepgRight[4];
     unsigned int i1, i2, i3, i4, i5, i6, i7;
     for( i1 = 1; i1 != m_iNumber-5; i1++)
         for(i2 = i1+1; i2 != m_iNumber-4; i2++)
             for(i3 = i2+1; i3 != m_iNumber-3; i3++)
                 for(i4 = i3+1; i4 != m_iNumber-2; i4++)
                     for(i5 = i4+1; i5 != m_iNumber-1; i5++)
                         for(i6 = i5+1; i6 != m_iNumber; i6++)
                             for(i7 = i6+1; i7 != m_iNumber+1; i7++)
     {         		   
		   		   		   
		   // Split of type 1: 321|7654
           SepgLeft[0] = i3;
           SepgLeft[1] = i2;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i5;
           SepgRight[3] = i4;
           dWei = MinSeptet(SepgLeft, SepgRight);
            
           CSeptet SepTemp1(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp1.GetIndex()] = SepTemp1;
           
           // Split of type 2: 421|7653
           SepgLeft[0] = i4;
           SepgLeft[1] = i2;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i5;
           SepgRight[3] = i3;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp2(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp2.GetIndex()] = SepTemp2;
           
           // Split of type 3: 431|7652
           SepgLeft[0] = i4;
           SepgLeft[1] = i3;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i5;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp3(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp3.GetIndex()] = SepTemp3;
           
           // Split of type 4: 432|7651
           SepgLeft[0] = i4;
           SepgLeft[1] = i3;
           SepgLeft[2] = i2;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i5;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp4(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp4.GetIndex()] = SepTemp4;
           
           // Split of type 5: 521|7643
           SepgLeft[0] = i5;
           SepgLeft[1] = i2;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i4;
           SepgRight[3] = i3;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp5(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp5.GetIndex()] = SepTemp5;
           
           // Split of type 6: 531|7642
           SepgLeft[0] = i5;
           SepgLeft[1] = i3;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i4;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp6(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp6.GetIndex()] = SepTemp6;
           
           // Split of type 7: 532|7641
           SepgLeft[0] = i5;
           SepgLeft[1] = i3;
           SepgLeft[2] = i2;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i4;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp7(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp7.GetIndex()] = SepTemp7;
           
           
           // Split of type 8: 541|7632
           SepgLeft[0] = i5;
           SepgLeft[1] = i4;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i3;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp8(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp8.GetIndex()] = SepTemp8;
           
           // Split of type 9: 542|7631
           SepgLeft[0] = i5;
           SepgLeft[1] = i4;
           SepgLeft[2] = i2;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i3;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp9(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp9.GetIndex()] = SepTemp9;
           
           // Split of type 10: 543|7621
           SepgLeft[0] = i5;
           SepgLeft[1] = i4;
           SepgLeft[2] = i3;
           SepgRight[0] = i7;
           SepgRight[1] = i6;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp10(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp10.GetIndex()] = SepTemp10;
           
           //Split of type 11: 621|7543
           SepgLeft[0] = i6;
           SepgLeft[1] = i2;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i5;
           SepgRight[2] = i4;
           SepgRight[3] = i3;
           dWei = MinSeptet(SepgLeft, SepgRight);
           CSeptet SepTemp11(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp11.GetIndex()] = SepTemp11;
           
           // Split of type 12: 631|7542
           SepgLeft[0] = i6;
           SepgLeft[1] = i3;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i5;
           SepgRight[2] = i4;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp12(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp12.GetIndex()] = SepTemp12;
           
           // Split of type 13: 632|7541
           SepgLeft[0] = i6;
           SepgLeft[1] = i3;
           SepgLeft[2] = i2;
           SepgRight[0] = i7;
           SepgRight[1] = i5;
           SepgRight[2] = i4;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp13(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp13.GetIndex()] = SepTemp13;
           
           // Split of type 14: 641|7532
           SepgLeft[0] = i6;
           SepgLeft[1] = i4;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i5;
           SepgRight[2] = i3;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp14(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp14.GetIndex()] = SepTemp14;
           
           // Split of type 15: 642|7531
           SepgLeft[0] = i6;
           SepgLeft[1] = i4;
           SepgLeft[2] = i2;
           SepgRight[0] = i7;
           SepgRight[1] = i5;
           SepgRight[2] = i3;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp15(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp15.GetIndex()] = SepTemp15;
                            
           // Split of type 16: 643|7521
           SepgLeft[0] = i6;
           SepgLeft[1] = i4;
           SepgLeft[2] = i3;
           SepgRight[0] = i7;
           SepgRight[1] = i5;
           SepgRight[2] = i2;
           SepgRight[3] = i1;           		
		   
		   dWei = MinSeptet(SepgLeft, SepgRight);		   
            
           CSeptet SepTemp16(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp16.GetIndex()] = SepTemp16;
           
           // Split of type 19: 653|7421
           SepgLeft[0] = i6;
           SepgLeft[1] = i5;
           SepgLeft[2] = i3;
           SepgRight[0] = i7;
           SepgRight[1] = i4;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp19(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp19.GetIndex()] = SepTemp19;
           
           // Split of type 18: 652|7431
           SepgLeft[0] = i6;
           SepgLeft[1] = i5;
           SepgLeft[2] = i2;
           SepgRight[0] = i7;
           SepgRight[1] = i4;
           SepgRight[2] = i3;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp18(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp18.GetIndex()] = SepTemp18;
           
           // Split of type 17: 651|7432
           SepgLeft[0] = i6;
           SepgLeft[1] = i5;
           SepgLeft[2] = i1;
           SepgRight[0] = i7;
           SepgRight[1] = i4;
           SepgRight[2] = i3;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp17(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp17.GetIndex()] = SepTemp17;
           
           // Split of type 20: 654|7321
           SepgLeft[0] = i6;
           SepgLeft[1] = i5;
           SepgLeft[2] = i4;
           SepgRight[0] = i7;
           SepgRight[1] = i3;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp20(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp20.GetIndex()] = SepTemp20;
           
           //Split of type 21: 721|6543
           SepgLeft[0] = i7;
           SepgLeft[1] = i2;
           SepgLeft[2] = i1;
           SepgRight[0] = i6;
           SepgRight[1] = i5;
           SepgRight[2] = i4;
           SepgRight[3] = i3;
           dWei = MinSeptet(SepgLeft, SepgRight);
           CSeptet SepTemp21(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp21.GetIndex()] = SepTemp21;
           
           // Split of type 22: 731|6542
           SepgLeft[0] = i7;
           SepgLeft[1] = i3;
           SepgLeft[2] = i1;
           SepgRight[0] = i6;
           SepgRight[1] = i5;
           SepgRight[2] = i4;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp22(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp22.GetIndex()] = SepTemp22;
           
           // Split of type 23: 732|6541
           SepgLeft[0] = i7;
           SepgLeft[1] = i3;
           SepgLeft[2] = i2;
           SepgRight[0] = i6;
           SepgRight[1] = i5;
           SepgRight[2] = i4;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp23(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp23.GetIndex()] = SepTemp23;
           
           // Split of type 25: 742|6531
           SepgLeft[0] = i7;
           SepgLeft[1] = i4;
           SepgLeft[2] = i2;
           SepgRight[0] = i6;
           SepgRight[1] = i5;
           SepgRight[2] = i3;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp25(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp25.GetIndex()] = SepTemp25;
           
           // Split of type 24: 741|6532
           SepgLeft[0] = i7;
           SepgLeft[1] = i4;
           SepgLeft[2] = i1;
           SepgRight[0] = i6;
           SepgRight[1] = i5;
           SepgRight[2] = i3;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp24(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp24.GetIndex()] = SepTemp24;
           
           // Split of type 26: 743|6521
           SepgLeft[0] = i7;
           SepgLeft[1] = i4;
           SepgLeft[2] = i3;
           SepgRight[0] = i6;
           SepgRight[1] = i5;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp26(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp26.GetIndex()] = SepTemp26;
           
           //Split of type 27: 751|6432
           SepgLeft[0] = i7;
           SepgLeft[1] = i5;
           SepgLeft[2] = i1;
           SepgRight[0] = i6;
           SepgRight[1] = i4;
           SepgRight[2] = i3;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);
            
           CSeptet SepTemp27(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp27.GetIndex()] = SepTemp27;
           
           // Split of type 28: 752|6431
           SepgLeft[0] = i7;
           SepgLeft[1] = i5;
           SepgLeft[2] = i2;
           SepgRight[0] = i6;
           SepgRight[1] = i4;
           SepgRight[2] = i3;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp28(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp28.GetIndex()] = SepTemp28;
           
           // Split of type 29: 753|6421
           SepgLeft[0] = i7;
           SepgLeft[1] = i5;
           SepgLeft[2] = i3;
           SepgRight[0] = i6;
           SepgRight[1] = i4;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp29(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp29.GetIndex()] = SepTemp29;
           
           // Split of type 30: 754|6321
           SepgLeft[0] = i7;
           SepgLeft[1] = i5;
           SepgLeft[2] = i4;
           SepgRight[0] = i6;
           SepgRight[1] = i3;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp30(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp30.GetIndex()] = SepTemp30;
           
           // Split of type 34: 764|5321
           SepgLeft[0] = i7;
           SepgLeft[1] = i6;
           SepgLeft[2] = i4;
           SepgRight[0] = i5;
           SepgRight[1] = i3;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp34(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp34.GetIndex()] = SepTemp34;
           
           // Split of type 33: 763|5421
           SepgLeft[0] = i7;
           SepgLeft[1] = i6;
           SepgLeft[2] = i3;
           SepgRight[0] = i5;
           SepgRight[1] = i4;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp33(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp33.GetIndex()] = SepTemp33;
           
           // Split of type 32: 762|5431
           SepgLeft[0] = i7;
           SepgLeft[1] = i6;
           SepgLeft[2] = i2;
           SepgRight[0] = i5;
           SepgRight[1] = i4;
           SepgRight[2] = i3;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp32(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp32.GetIndex()] = SepTemp32;
           
           // Split of type 31: 761|5432
           SepgLeft[0] = i7;
           SepgLeft[1] = i6;
           SepgLeft[2] = i1;
           SepgRight[0] = i5;
           SepgRight[1] = i4;
           SepgRight[2] = i3;
           SepgRight[3] = i2;
           dWei = MinSeptet(SepgLeft, SepgRight);   
           CSeptet SepTemp31(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp31.GetIndex()] = SepTemp31;
           
           //Split of type 35: 765|4321 
           SepgLeft[0] = i7;
           SepgLeft[1] = i6;
           SepgLeft[2] = i5;
           SepgRight[0] = i4;
           SepgRight[1] = i3;
           SepgRight[2] = i2;
           SepgRight[3] = i1;
           dWei = MinSeptet(SepgLeft, SepgRight);  
           CSeptet SepTemp35(SepgLeft, SepgRight, dWei);
           m_gSeptet[SepTemp35.GetIndex()] = SepTemp35;
     }
 }
