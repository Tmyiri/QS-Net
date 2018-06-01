#include <iostream>
#include "CSplitSys.h"

using namespace std;

// calculate w(21|543) 

// Method 1 w(21|543)=  min{
//                      w(43|21)-w(51|43)+w(51|32)-w(54|32)+w(54|21),  
//                      w(53|21)-w(53|41)+w(41|32)-w(54|32)+w(54|21), 
//                      w(43|21)-w(51|43)+w(51|42)-w(53|42)+w(53|21),   
//                      w(54|21)-w(54|31)+w(42|31)-w(53|42)+w(53|21),  
//                      w(53|21)-w(53|41)+w(52|41)-w(52|43)+w(43|21),           
//                      w(54|21)-w(54|31)+w(52|31)-w(52|43)+w(43|21).}

 double CSplitSys::MinQuintet(unsigned int gLeft[2], unsigned int gRight[3])
{  
   CQuartet qTemp;
   
   // specify all the possible blocks
   unsigned int iB21[2];
   iB21[0] = gLeft[0]; 
   iB21[1] = gLeft[1];
   
   unsigned int  iB31[2];
   iB31[0] = gRight[2]; 
   iB31[1] = gLeft[1];
   
   unsigned int iB32[2];
   iB32[0] = gRight[2]; 
   iB32[1] = gLeft[0];
   
   unsigned int iB41[2];
   iB41[0] = gRight[1]; 
   iB41[1] = gLeft[1];
   
   unsigned int iB42[2];
   iB42[0] = gRight[1]; 
   iB42[1] = gLeft[0];
   
   unsigned int iB43[2];
   iB43[0] = gRight[1];
   iB43[1] = gRight[2];
   
   unsigned int iB51[2];
   iB51[0] = gRight[0];
   iB51[1] = gLeft[1];
   
   unsigned int iB52[2];
   iB52[0] = gRight[0];
   iB52[1] = gLeft[0];
   
   unsigned int iB53[2];
   iB53[0] = gRight[0];
   iB53[1] = gRight[2];
   
   unsigned int iB54[2];
   iB54[0] = gRight[0];
   iB54[1] = gRight[1];
   
  // First: w(43|21)-w(51|43)+w(51|32)-w(54|32)+w(54|21)
    
   double dTemp = 0;
   // + w(43|21)
   qTemp.SetLabel(iB43, iB21);
   double w4321 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w4321;
   
   //- w(51|43)
   qTemp.SetLabel(iB51, iB43);
   double w5143 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp -= w5143;
   
   //+w(51|32)
   qTemp.SetLabel(iB51, iB32);
   double w5132 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w5132;
   
   //-w(54|32)
   qTemp.SetLabel(iB54, iB32);
   double w5432 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp -= w5432;
   
   //+w(54|21)
   qTemp.SetLabel(iB54, iB21);
   double w5421 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w5421;
   
   if(dTemp <= 0) return 0;
   
   double dWei =dTemp;
   
   // Second: w(53|21)-w(53|41)+w(41|32)-w(54|32)+w(54|21)
   dTemp=0;
   // + w(53|21)
   qTemp.SetLabel(iB53, iB21);
   double w5321 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w5321;
   
   //- w(53|41)
   qTemp.SetLabel(iB53, iB41);
   double w5341 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp -= w5341;
   
   //+w(41|32)
   qTemp.SetLabel(iB41, iB32);
   double w4132 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w4132;
   
   //-w(54|32)
   dTemp -= w5432;
   
   //+w(54|21)
   dTemp += w5421;
   
   if(dTemp <= 0) return 0;
   if(dWei > dTemp) dWei = dTemp;
   
   // Third: w(43|21)-w(51|43)+w(51|42)-w(53|42)+w(53|21)
   dTemp=0;
   // + w(43|21)
   dTemp += w4321;
   
   //- w(51|43)
   dTemp -= w5143;
   
   //+w(51|42)
   qTemp.SetLabel(iB51, iB42);
   double w5142 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w5142;
   
   //-w(53|42)
   qTemp.SetLabel(iB53, iB42);
   double w5342 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp -= w5342;
   
   //+w(53|21)
   dTemp += w5321;
   
   if(dTemp <= 0) return 0;
   if(dWei > dTemp) dWei = dTemp;
    
   // Fourth: w(54|21)-w(54|31)+w(42|31)-w(53|42)+w(53|21),
   dTemp=0;
   // + w(54|21)
   dTemp += w5421;
   
   //-w(54|31)
   qTemp.SetLabel(iB54, iB31);
   double w5431 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp -= w5431;
   
   //+w(42|31)
   qTemp.SetLabel(iB42, iB31);
   double w4231 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w4231;
   
   //-w(53|42)
   dTemp -= w5342;
   
   //+w(53|21)
   dTemp += w5321;
   
   if(dTemp <= 0) return 0;
   if(dWei > dTemp) dWei = dTemp;
   
   // Fifth: w(53|21)-w(53|41)+w(52|41)-w(52|43)+w(43|21), 
   dTemp=0;
   // + w(53|21)
   dTemp += w5321;
   
   //- w(53|41)
   dTemp -= w5341;
   
   //+w(52|41)
   qTemp.SetLabel(iB52, iB41);
   double w5241 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w5241;
   
   //-w(52|43)
   qTemp.SetLabel(iB52, iB43);
   double w5243 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp -= w5243;
   
   //+w(43|21)
   dTemp += w4321;
   
   if(dTemp <= 0) return 0;
   if(dWei > dTemp) dWei = dTemp;
   
   // Sixth: w(54|21)-w(54|31)+w(52|31)-w(52|43)+w(43|21)
   dTemp=0;
   // + w(54|21)
   dTemp += w5421;
   
   //- w(54|31)
   dTemp -= w5431;
   
   //+w(52|31)
   qTemp.SetLabel(iB52, iB31);
   double w5231 = m_gQuartet[qTemp.GetIndex()].GetWeight();
   dTemp += w5231;
   
   //-w(52|43)
   dTemp -= w5243;
   
   //+w(43|21)
   dTemp += w4321;
   
   if(dTemp <= 0) return 0;
   if(dWei > dTemp) dWei = dTemp;
   
   return 0.5*dWei;
}

// EOF
