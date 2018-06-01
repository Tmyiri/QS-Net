#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

#include "CSplitSys.h"
#include "auxilary.h"

using namespace std;

// Save data from multiple sequence alignment file: fasta format
void CSplitSys::ReadAlignment(const string& InputFile)
{
     // Open file filename for reading
     ifstream fin(InputFile.c_str());                                 
               
     if(!fin ) 
     {   
         cerr << "Error to  open " << InputFile.c_str() << " for reading!" << endl;  
         system("pause"); 
         exit(1);  
     }
     else
     {    
        string sLine;                   // sLine stores the strings in a line
		vector<string> gSeq;            // stores the sequences of the taxa set
				
        string sSeq = "";
                
	     // Read the file line by line 
         while(getline(fin, sLine) )                                         
         {   
             if(sLine == "") continue;         // If it is a blank line, go back and read the next line
             else
		     {
                 if(sLine[0] == '>')           // Indicate that the sequence of a new species comes 
                 { 
                      gSeq.push_back(sSeq);    // push back the sequence of the previous species to gSeq
                 
					  m_gName.push_back(sLine.substr(1, sLine.size()-1));
						  
                      sSeq = "";
                 }    
                 else
                 {
				     sLine = trim(sLine);      // remove the possible spaces at the begining and end of a line
                     sSeq += sLine;
                 }                
             }
         }
		
        gSeq.push_back(sSeq);                  // push the last sequence
        m_iNumber = m_gName.size();
		
		
	    // Calculate triplet weights from the multiple alignment 

		// Set all the 3*(n,3)  triplets to be empty triplets;
 	    CTriplet tEmpty;
        for(unsigned long long n = 0; n != 3*Binom(m_iNumber, 3); ++n)
             m_gTriplet.push_back(tEmpty);
        
        unsigned int gRight[2];    
        
		string::size_type i, j, k, l, i1, i2, i3, i4, i5, i6; 
		
		for(i = 1; i != m_iNumber-1; ++i)
		   for( j = i+1; j != m_iNumber; ++j)
			   for( k = j+1; k != m_iNumber +1; ++k)
		{
		      // initialize 3 triplets 
			  
			  // Triplet 1: i|kj, e.g. 1|32
			  gRight[0] = k;
			  gRight[1] = j;
			  
		      double dWei = (double) CalTriplet(gSeq[i], gSeq[k], gSeq[j]);
		      CTriplet tOne(i, gRight, dWei);
		      
		      // push triplet 1 to the corresponding position of  m_gTriplet
			  m_gTriplet[tOne.GetIndex()] = tOne;
			  
			  // push triplet 2: j|ki to m_gTriplet 2|31
			  dWei = (double) CalTriplet(gSeq[j], gSeq[k], gSeq[i]);
              gRight[1] = i;
              
		      CTriplet tTwo(j, gRight, dWei);
			  m_gTriplet[tTwo.GetIndex()] = tTwo;
			  
			  // push triplet 3: k|ji to m_gTriplet 3|21
			  dWei = (double) CalTriplet(gSeq[k], gSeq[j], gSeq[i]);
              gRight[0] = j;
              
		      CTriplet tThree(k, gRight, dWei);
			  m_gTriplet[tThree.GetIndex()] = tThree;
		}
		
	    // Calculate quartet weights from the multiple alignment 	
        
	    // specify all 3(n,4) quartets to be empty quartets
       	CQuartet qEmpty;
        for(unsigned long long n = 0; n != 3*Binom(m_iNumber,4); ++n)
             m_gQuartet.push_back(qEmpty);
       
        unsigned int gLeft[2];      
       
        for(i = 1; i != m_iNumber-2; ++i)
		   for( j = i+1;j != m_iNumber-1; ++j)
			  for( k = j+1; k != m_iNumber; ++k)
                 for( l = k+1; l != m_iNumber+1; ++l)
        {
                 // Spcify quartet 1: li|kj and its weight w(li|kj)
                 gLeft[0] = l;
                 gLeft[1] = i;
                 gRight[0] = k;
                 gRight[1] = j;
                 double dWei = (double) CalQuartet(gSeq[l], gSeq[i], gSeq[k], gSeq[j]);
                 CQuartet qOne(gLeft, gRight, dWei);
                 m_gQuartet[qOne.GetIndex()] = qOne;
                         
                 // Spcify quartet 2: lj|ki and its weight w(lj|ki)
                 swap(gLeft[1], gRight[1]);
                 dWei = (double)CalQuartet(gSeq[l], gSeq[j], gSeq[k], gSeq[i]);
                 CQuartet qTwo(gLeft, gRight, dWei);
                 m_gQuartet[qTwo.GetIndex()] = qTwo;

                 // Spcify quartet 3: lk|ji and its weight w(lk|ji)
                 swap(gLeft[1], gRight[0]);
                 dWei=(double)CalQuartet(gSeq[l], gSeq[k], gSeq[j], gSeq[i]);
                 CQuartet qThree(gLeft, gRight, dWei);
                 m_gQuartet[qThree.GetIndex()] = qThree;
        } 
          
		// Calculate sextet weights from the multiple sequence alignment 	
        
	    // specify all 10(n,6) quartets to be empty sextets       
        CSextet sEmpty;
        for(unsigned long long n = 0; n != 10 * Binom(m_iNumber, 6); ++n)
             m_gSextet.push_back(sEmpty);
       
        unsigned int sLeft[3]; 
		unsigned int sRight[3];     
       
        for(i1 = 1; i1 != m_iNumber-4; ++i1)
		   for( i2 = i1+1;i2 != m_iNumber-3; ++i2)
			 for( i3 = i2+1; i3 != m_iNumber-2; ++i3)
                 for( i4 = i3+1; i4 != m_iNumber-1; ++i4)
                    for(i5 = i4+1; i5 != m_iNumber; ++i5)
                       for(i6 = i5+1; i6 != m_iNumber+1; ++i6)
         {
	           // e.g 1: 621|543
	           sLeft[0] = i6;
	           sLeft[1] = i2;
	           sLeft[2] = i1;
	           sRight[0] = i5;
	           sRight[1] = i4;
	           sRight[2] = i3;
	           double dWei = (double) CalSextet(gSeq[i6], gSeq[i2], gSeq[i1], gSeq[i5], gSeq[i4], gSeq[i3]);
               CSextet sOne(sLeft, sRight, dWei);
               m_gSextet[sOne.GetIndex()] = sOne;
	           
	           // e.g 2: 631|542
	           sLeft[0] = i6;
	           sLeft[1] = i3;
	           sLeft[2] = i1;
	           sRight[0] = i5;
	           sRight[1] = i4;
	           sRight[2] = i2;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i3], gSeq[i1], gSeq[i5], gSeq[i4], gSeq[i2]);
               CSextet sTwo(sLeft, sRight, dWei);
               m_gSextet[sTwo.GetIndex()] = sTwo;
	           
	           // e.g 3: 632|541
	           sLeft[0] = i6;
	           sLeft[1] = i3;
	           sLeft[2] = i2;
	           sRight[0] = i5;
	           sRight[1] = i4;
	           sRight[2] = i1;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i3], gSeq[i2], gSeq[i5], gSeq[i4], gSeq[i1]);
               CSextet sThree(sLeft, sRight, dWei);
               m_gSextet[sThree.GetIndex()] = sThree;
	           
	           // e.g 5: 642|531
	           sLeft[0] = i6;
	           sLeft[1] = i4;
	           sLeft[2] = i2;
	           sRight[0] = i5;
	           sRight[1] = i3;
	           sRight[2] = i1;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i4], gSeq[i2], gSeq[i5], gSeq[i3], gSeq[i1]);
               CSextet sFive(sLeft, sRight, dWei);
               m_gSextet[sFive.GetIndex()] = sFive;
	           
	           // e.g 4: 641|532
	           sLeft[0] = i6;
	           sLeft[1] = i4;
	           sLeft[2] = i1;
	           sRight[0] = i5;
	           sRight[1] = i3;
	           sRight[2] = i2;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i4], gSeq[i1], gSeq[i5], gSeq[i3], gSeq[i2]);
               CSextet sFour(sLeft, sRight, dWei);
               m_gSextet[sFour.GetIndex()] = sFour;
	           
	           // e.g 6: 643|521
	           sLeft[0] = i6;
	           sLeft[1] = i4;
	           sLeft[2] = i3;
	           sRight[0] = i5;
	           sRight[1] = i2;
	           sRight[2] = i1;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i4], gSeq[i3], gSeq[i5], gSeq[i2], gSeq[i1]);
               CSextet sSix(sLeft, sRight, dWei);
               m_gSextet[sSix.GetIndex()] = sSix;
	           
	           // e.g 9: 653|421
	           sLeft[0] = i6;
	           sLeft[1] = i5;
	           sLeft[2] = i3;
	           sRight[0] = i4;
	           sRight[1] = i2;
	           sRight[2] = i1;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i5], gSeq[i3], gSeq[i4], gSeq[i2], gSeq[i1]);
               CSextet sNine(sLeft, sRight, dWei);
               m_gSextet[sNine.GetIndex()] = sNine;
	           
	           // e.g 8: 652|431
	           sLeft[0] = i6;
	           sLeft[1] = i5;
	           sLeft[2] = i2;
	           sRight[0] = i4;
	           sRight[1] = i3;
	           sRight[2] = i1;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i5], gSeq[i2], gSeq[i4], gSeq[i3], gSeq[i1]);
               CSextet sEight(sLeft, sRight, dWei);
               m_gSextet[sEight.GetIndex()] = sEight;
	           
	           // e.g 7: 651|432
	           sLeft[0] = i6;
	           sLeft[1] = i5;
	           sLeft[2] = i1;
	           sRight[0] = i4;
	           sRight[1] = i3;
	           sRight[2] = i2;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i5], gSeq[i1], gSeq[i4], gSeq[i3], gSeq[i2]);
               CSextet sSeven(sLeft, sRight, dWei);
               m_gSextet[sSeven.GetIndex()] = sSeven;
	           
	          // e.g 10: 654|321
	           sLeft[0] = i6;
			   sLeft[1] = i5;
			   sLeft[2] = i4;
	           sRight[0] = i3;
	           sRight[1] = i2;
	           sRight[2] = i1;
	           dWei = (double) CalSextet(gSeq[i6], gSeq[i5], gSeq[i4], gSeq[i3], gSeq[i2], gSeq[i1]);
               CSextet sTen(sLeft, sRight, dWei);
               m_gSextet[sTen.GetIndex()] = sTen;
        } 
	}
}
  
// EOF
