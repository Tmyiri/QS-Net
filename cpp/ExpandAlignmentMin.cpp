#include <string>
#include <iostream>
#include <fstream>
#include "CSplitSys.h"
#include "auxilary.h"

using namespace std;

// Expand and generate the complete splits from the fasta format file  
void CSplitSys::ExpandAlignmentMin()
{
   // Define gActive to store active splits and gTemp to store temp splits           
   vector<CSplit> gActive, gTemp;    
   vector<unsigned int> gLeft, gRight, gCount, gCount2;
   unsigned int i, k; 
   unsigned int it1, it2; 
   double dWei;           
   
   vector<unsigned int> gLt(2,0),gRf(3,0);
   
   //take 1£¬2£¬....£¬ m_iNumber into gCount
   for(i = m_iNumber; i != 0; i--){
   	  gCount.push_back(i);
   }
   
   //generate 21|m_iNumber m_iNumber-1..3, 31| m_iNumber m_iNumber-1....2, ... ,
   for(it1 = 1; it1 != m_iNumber; it1++ ) 
      for(it2 = it1 + 1; it2 != m_iNumber+1; it2++) 
      {
      	 vector<unsigned int> gCount1;
      	 gLt[0] = it2;
      	 gLt[1] = it1;
      	 
      	 gCount1 = gCount;
      	 vcomplete(gCount1, it1, it2);
      	 
      	 dWei = MinWei2n( gLt, gCount1);
      	 CSplit sOne(gLt, gCount1, dWei);
         m_gSplit.push_back(sOne);
	  }
	 
  
   //take 1,2,....,7 into gCount2 
   for(i = 7; i != 0; i--){
   	  gCount2.push_back(i);
   }
    
   unsigned int gLt1[3], gRf1[5];
   int i1; 
   //generate 821|76..3  831|76..2 ...
   for(it1 = 1; it1 != 7; it1++ ) 
      for(it2 = it1 + 1; it2 != 8; it2++) 
      {
      	 vector<unsigned int> gCount3;
      	 gRf[0] = 8;
      	 gRf[1] = it2;
      	 gRf[2] = it1;
      	 
      	 gLt1[0] = gRf[0];
      	 gLt1[1] = gRf[1];
      	 gLt1[2] = gRf[2];
      	 
      	 gCount3 = gCount2;
      	 vcomplete(gCount3, it1, it2);
      	 for(i1 = 0; i1 < 5; i1++){
      	 	gRf1[i1] = gCount3[i1];
		 }
      	 
      	 dWei = MinWei35( gLt1, gRf1);
      	 CSplit sTemp(gRf, gCount3, dWei);
         gActive.push_back(sTemp);
	  }
     
	
	//take 8 to the left of 321|7654 421|7653 431|7652... or the right 
	unsigned int gLL[4];
	vector<CSeptet>::iterator it34;
     for( it34 = m_gSeptet.begin(); it34 != m_gSeptet.begin()+35; ++it34)
   { 	 
         vector<unsigned int> gLone, gRone;
         unsigned int *gL = (*it34).GetLblk();
         unsigned int *gR = (*it34).GetRblk();
         
		// add 8 to the left
		gLL[0] = 8;
		gLone.push_back(8);
		gLL[1] = gL[0];
		gLone.push_back(gL[0]);
		gLL[2] = gL[1];
		gLone.push_back(gL[1]);
		gLL[3] = gL[2];
		gLone.push_back(gL[2]);
		
		gRone.push_back(gR[0]);
		gRone.push_back(gR[1]);
		gRone.push_back(gR[2]);
		gRone.push_back(gR[3]);

		dWei = MinWei44(gLL, gR);
		CSplit sLeft(gLone, gRone, dWei);
		gActive.push_back(sLeft);	
		
		// add 8 to the right 
		vector<unsigned int> gLtwo;
		gLtwo.push_back(gL[0]);
		gLtwo.push_back(gL[1]);
		gLtwo.push_back(gL[2]);
		
		gRf1[0] = 8;
		gRf1[1] = gR[0];
		gRf1[2] = gR[1];
		gRf1[3] = gR[2];
		gRf1[4] = gR[3];
		
		gRone.insert(gRone.begin(),8);
		
		dWei = MinWei35(gL, gRf1);
		CSplit sRight(gLtwo, gRone, dWei);
		gActive.push_back(sRight);
	}
    	
    // delete the CSplits with weight 0
	ShrinkSplits(gActive); 
	
	vector<CSplit>::iterator iter;
    // Case: k >= 9
	for(k = 9; k != m_iNumber+1; ++k)
	{
        cout << "Step " << k << " Starts!" << endl;
        
		// As the same, there are two steps:
		// Step 1: generating type k1|(k-1)...2, k2|..., k3|..., k(k-1)|... and their weights 
		
		
		//take k-1 into gCount2
        gCount2.insert(gCount2.begin(), k-1);

       //generate 921|876..3  931|76..2 ...
       for(it1 = 1; it1 != k-1; it1++ ) 
	      for(it2 = it1 + 1; it2 != k; it2++) 
	      {
	      	 vector<unsigned int> gLb;
	      	 vector<unsigned int> gCount4;
	      	 gLb.push_back(k);
	      	 gLb.push_back(it2);
	      	 gLb.push_back(it1);
	      	 
	      	 gCount4 = gCount2;
	      	 vcomplete(gCount4, it1, it2);
	      	 
	      	 dWei = MinWei3n( gLb, gCount4);
	      	 CSplit sFirst(gLb, gCount4, dWei);
	         gTemp.push_back(sFirst);
		  }
		 
		 // Step 2: generating other types by adding k to the left and right hand side of active_CSplits 
		 
		 for( iter = gActive.begin(); iter != gActive.end(); ++iter)
    	{          
             // Add k to the left
	    	gLeft = (*iter).GetLblk();
	    	gRight = (*iter).GetRblk();
	    	gLeft.insert(gLeft.begin(), k);
		            
	    	dWei = MinWeimn(gLeft, gRight);
                 
	     	CSplit sLeft1(gLeft, gRight, dWei);
		    gTemp.push_back(sLeft1);
		
		   // Add k to the right
	    	gLeft = (*iter).GetLblk();
	    	gRight = (*iter).GetRblk();
	    	gRight.insert(gRight.begin(), k);
		
            if(gLeft.size() == 3)
	             dWei = MinWei3n(gLeft, gRight);
            else
                 dWei = MinWeimn(gLeft, gRight);
                 
	     	CSplit sRight1(gLeft, gRight, dWei);
		    gTemp.push_back(sRight1);
	    }
	    
        // Delete the CSplits with weight 0
	    ShrinkSplits(gTemp);
	    gActive = gTemp;
	    gTemp.clear();
	    
	    cout<<"Step "<<k<<" Ended!"<<endl;
    }
	for( iter = gActive.begin(); iter != gActive.end(); ++iter)
	    m_gSplit.push_back(*iter);
	ShrinkSplits(m_gSplit);
}
