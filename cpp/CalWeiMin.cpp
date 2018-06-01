#include "CSplitSys.h"
#include "auxilary.h"
#include "CQuintet.h"
#include "CSeptet.h"
#include<algorithm> 


// Calculate w(321|87654)
double CSplitSys::MinWei35(unsigned int gLeft[3], unsigned int gRight[5])
{  
   CSeptet tTemp;
  // 321|8765
   unsigned int gRb[4];
   gRb[0] = gRight[0];
   gRb[1] = gRight[1];
   gRb[2] = gRight[2];
   gRb[3] = gRight[3];
   tTemp.SetLabel(gLeft, gRb);
   double dMin = m_gSeptet[tTemp.GetIndex()].GetWeight();
   
   if(dMin <= 0) return 0;
   
   // 321|8764
   gRb[3] = gRight[4];
   tTemp.SetLabel(gLeft, gRb);
   double dWei = m_gSeptet[tTemp.GetIndex()].GetWeight();
			 
   if(dWei <= 0) return 0;
   else if(dMin > dWei) dMin = dWei;
   
   // 321|8754
   gRb[2] = gRight[3];
   tTemp.SetLabel(gLeft, gRb);
   dWei = m_gSeptet[tTemp.GetIndex()].GetWeight();
			 
   if(dWei <= 0) return 0;
   else if(dMin > dWei) dMin = dWei;
   
   // 321|8654
   gRb[1] = gRight[2];
   tTemp.SetLabel(gLeft, gRb);
   dWei = m_gSeptet[tTemp.GetIndex()].GetWeight();
			 
   if(dWei <= 0) return 0;
   else if(dMin > dWei) dMin = dWei;
   
   // 321|7654
   gRb[0] = gRight[1];
   tTemp.SetLabel(gLeft, gRb);
   dWei = m_gSeptet[tTemp.GetIndex()].GetWeight();
			 
   if(dWei <= 0) return 0;
   else if(dMin > dWei) dMin = dWei;
         
   return dMin;
}

// Calculate w(21|n...6543)
double CSplitSys::MinWei2n(const vector<unsigned int>& gLeft,const vector<unsigned int>& gRight)
{  
   double dMin = 10000000;
   double dWei;
   CQuintet tTemp;
   unsigned int gLb[2], gRb[3];
   gLb[0] = gLeft[0];
   gLb[1] = gLeft[1];
   
   unsigned int it1, it2, it3;

   for( it1 = 0; it1 != gRight.size()-2; ++it1)
      for(it2 = it1+1; it2 != gRight.size()-1; ++it2)
	     for(it3 = it2+1; it3 != gRight.size(); ++it3)
		 {
			 gRb[0] = gRight[it1];
			 gRb[1] = gRight[it2];
			 gRb[2] = gRight[it3];
			 tTemp.SetLabel(gLb, gRb);
			
			 dWei = m_gQuintet[tTemp.GetIndex()].GetWeight();
			 
             if(dWei <= 0) return 0;
             else if(dMin > dWei) dMin = dWei;
         }
   return dMin;
}

//  calculate w(1234|5678) 
double CSplitSys::MinWei44(unsigned int gLeft[4], unsigned int gRight[4])
{
   double dWei, dTemp;
   unsigned int temp;
   CSeptet tTemp;
   
   // case1: w(123|5678)-w(123|45678)
   unsigned int gLb[3];
   gLb[0] = gLeft[0];
   gLb[1] = gLeft[1];
   gLb[2] = gLeft[2];
   unsigned int gRb[5];
   gRb[0] = gRight[0];
   gRb[1] = gRight[1];
   gRb[2] = gRight[2];
   gRb[3] = gRight[3];
   gRb[4] = gLeft[3];
   
   tTemp.SetLabel(gLb, gRight);
   dTemp = m_gSeptet[tTemp.GetIndex()].GetWeight();
   dTemp -= MinWei35( gLb, gRb);
   dWei = dTemp;
   if( dWei <= 0) return 0;
   
   //case 2: w(124|5678)-w(124|35678)
   swap(gLb[2], gRb[4]);
//   temp = gLb[2];
//   gLb[2] = gRb[4];
//   gRb[4] = temp;
   tTemp.SetLabel(gLb, gRight);
   dTemp = m_gSeptet[tTemp.GetIndex()].GetWeight();
   dTemp -= MinWei35( gLb, gRb);
   if(dTemp <= 0) return 0;
   else if(dWei > dTemp) dWei = dTemp;
   
   //case 3: w(134|5678)-w(124|25678)
   swap(gLb[1], gRb[4]);
//   temp = gLb[1];
//   gLb[1] = gRb[4];
//   gRb[4] = temp;
   tTemp.SetLabel(gLb, gRight);
   dTemp = m_gSeptet[tTemp.GetIndex()].GetWeight();
   dTemp -= MinWei35( gLb, gRb);
   if(dTemp <= 0) return 0;
   else if(dWei > dTemp) dWei = dTemp;
   
   //case 4: w(234|5678)-w(234|15678)
   swap(gLb[0], gRb[4]);
//   temp = gLb[0];
//   gLb[0] = gRb[4];
//   gRb[4] = temp;
   tTemp.SetLabel(gLb, gRight);
   dTemp = m_gSeptet[tTemp.GetIndex()].GetWeight();
   dTemp -= MinWei35( gLb, gRb);
   if(dTemp <= 0) return 0;
   else if(dWei > dTemp) dWei = dTemp;
        
   
   //case 5: w(567|1234)-w(567|12348)
   unsigned int gLnew[3];
   gLnew[0] = gRight[0];
   gLnew[1] = gRight[1];
   gLnew[2] = gRight[2];
   unsigned int gRnew[5];
   gRnew[0] = gLeft[0];
   gRnew[1] = gLeft[1];
   gRnew[2] = gLeft[2];
   gRnew[3] = gLeft[3];
   gRnew[4] = gRight[3];
   
   tTemp.SetLabel(gLnew, gLeft);
   dTemp = m_gSeptet[tTemp.GetIndex()].GetWeight();
   dTemp -= MinWei35( gLnew, gRnew);
   if(dTemp <= 0) return 0;
   else if(dWei > dTemp) dWei = dTemp;     
   
   //case 6: w(568|1234)-w(568|12347)
   swap(gLnew[2], gRnew[4]);
//   temp = gLnew[2];
//   gLnew[2] = gRnew[4];
//   gRnew[4] = temp;
   tTemp.SetLabel(gLnew, gLeft);
   dTemp = m_gSeptet[tTemp.GetIndex()].GetWeight();
   dTemp -= MinWei35( gLnew, gRnew);
   if(dTemp <= 0) return 0;
   else if(dWei > dTemp) dWei = dTemp; 
   
   //case 7: w(578|1234)-w(578|12346)
   swap(gLnew[1], gRnew[4]);
//   temp = gLnew[1];
//   gLnew[1] = gRnew[4];
//   gRnew[4] = temp;
   tTemp.SetLabel(gLnew, gLeft);
   dTemp = m_gSeptet[tTemp.GetIndex()].GetWeight();
   dTemp -= MinWei35( gLnew, gRnew);
   if(dTemp <= 0) return 0;
   else if(dWei > dTemp) dWei = dTemp;
   
   //case 8: w(678|1234)-w(678|12345)
   swap(gLnew[0], gRnew[4]);
//   temp = gLnew[0];
//   gLnew[0] = gRnew[4];
//   gRnew[4] = temp;
   tTemp.SetLabel(gLnew, gLeft);
   dTemp = m_gSeptet[tTemp.GetIndex()].GetWeight();
   dTemp -= MinWei35( gLnew, gRnew);
   if(dTemp <= 0) return 0;
   else if(dWei > dTemp) dWei = dTemp; 
   
   return dWei;
}

// calculate w(321|n...7654)
double CSplitSys::MinWei3n(const vector<unsigned int>& gLeft, const vector<unsigned int>& gRight)
{
   double dMin3n = 10000000;
   double dWei;
   CSeptet SepTemp;
   unsigned int gLb[3], gRb[4];
   gLb[0] = gLeft[0];
   gLb[1] = gLeft[1];
   gLb[2] = gLeft[2];
   
   unsigned int it1, it2, it3, it4;

   for( it1 = 0; it1 != gRight.size()-3; ++it1)
      for(it2 = it1+1; it2 != gRight.size()-2; ++it2)
	     for(it3 = it2+1; it3 != gRight.size()-1; ++it3)
	        for(it4 = it3+1; it4 != gRight.size(); ++it4)
		    {
				 gRb[0] = gRight[it1];
				 gRb[1] = gRight[it2];
				 gRb[2] = gRight[it3];
				 gRb[3] = gRight[it4];
				 SepTemp.SetLabel(gLb, gRb);
				 dWei = m_gSeptet[SepTemp.GetIndex()].GetWeight();
				 
	             if(dWei <= 0) return 0;
	             else if(dMin3n > dWei) dMin3n = dWei;
           }
   return dMin3n;
}


//  calculate w(m|n)
double CSplitSys::MinWeimn(const vector<unsigned int>& gLeft, const vector<unsigned int>& gRight)
{
   double dWei = 10000000;
   double dTemp; 
   unsigned int gLb[4], gRb[4];
   
   unsigned int it1, it2, it3, it4, it5, it6, it7 ,it8;
   
     for(it1 = 0; it1 != gLeft.size()-3; ++it1)
	   for(it2 = it1 + 1; it2 != gLeft.size()-2; ++it2)
	      for(it3 = it2 + 1; it3 != gLeft.size()-1; ++it3)
		     for(it4 = it3 + 1; it4 != gLeft.size(); ++it4)
			    for(it5 = 0; it5 != gRight.size()-3; ++it5)
   	  	 	        for(it6 = it5 + 1; it6 != gRight.size()-2; ++it6)
   	  	 		       for(it7 = it6 + 1; it7 != gRight.size()-1; ++it7)
   	  	 		          for(it8 = it7 + 1; it8 != gRight.size(); ++it8) 	  	 		     
		    {
				     	gLb[0] = gLeft[it1];
				    	gLb[1] = gLeft[it2];
					    gLb[2] = gLeft[it3];
				 	    gLb[3] = gLeft[it4];
					    gRb[0] = gRight[it5];	
					    gRb[1] = gRight[it6];
					    gRb[2] = gRight[it7];
					    gRb[3] = gRight[it8];					 
					
					    dTemp = MinWei44(gLb, gRb);
								 
					    if(dTemp <= 0) return 0;
					    else if(dWei > dTemp) dWei = dTemp;	 
	       } 
				        
   return dWei;
}

// EOF
