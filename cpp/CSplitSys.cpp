#include <string>
#include <vector> 
#include <algorithm>
#include <iostream>
#include <fstream>
#include "CSplitSys.h"
#include "auxilary.h"
                                             
// Operator constructor
CSplitSys& CSplitSys::operator=(const CSplitSys& rhs)
{
	m_iNumber = rhs.m_iNumber;
	m_gName = rhs.m_gName;
	m_gTriplet = rhs.m_gTriplet;
	m_gQuartet = rhs.m_gQuartet;
	m_gQuintet = rhs.m_gQuintet;
	m_gSextet = rhs.m_gSextet;
	m_gSeptet = rhs.m_gSeptet;
	m_gSplit = rhs.m_gSplit;
	return *this;
}   

// Function to order sextets
bool CmpSextet(const CSextet& sA, const CSextet& sB)
{
	return sA.GetIndex() < sB.GetIndex();
}

// Reorder the sextets
void CSplitSys::OrderSextets()
{
	sort(m_gSextet.begin(), m_gSextet.end(), CmpSextet);
}
// Function to order quartets
bool CmpInd(const CQuartet &qA, const CQuartet &qB) 
{ 
     return qA.GetIndex() < qB.GetIndex();
}

// Reorder the quartets
void CSplitSys::OrderQuartets()
{
     sort(m_gQuartet.begin(), m_gQuartet.end(), CmpInd);
 }                                                                                                                                     
 
// Function to order triplets
bool CmpTriplet(const CTriplet &tA, const CTriplet &tB) 
{ 
     return tA.GetIndex() < tB.GetIndex();
}

// Reorder the triplets
void CSplitSys::OrderTriplets()
{
     sort(m_gTriplet.begin(), m_gTriplet.end(), CmpTriplet);
 }           
 
 
// Delete the split in split_sys that has weight 0
void CSplitSys::ShrinkSplits(vector<CSplit>& sTemp)
{
    vector<CSplit>::iterator iter;
    for(iter = sTemp.begin(); iter != sTemp.end(); ++iter)
	   {
	      if((*iter).GetWeight() == 0) 
		  {
              sTemp.erase(iter);
		      --iter;
		  }
	   } 
}

//Output split with nonzero weight to file 
void CSplitSys::Output(string& OutputFile, double dPercent)
{
	// Output the final complete Splits and weights into the OutputFile
     ofstream fout;
	 
     fout.open( OutputFile.c_str(),   ios_base::out   |   ios_base::trunc);  
		    
     if(!fout) 
     {   
        cerr << "Error opening " << OutputFile << " for output!" << endl;  
        system("pause"); 
        exit(1);  
     }
     else
     {
	     // Output taxa number and labels
	     fout << "#nexus" << endl;
		 fout << endl;
		 
		 fout << "BEGIN Taxa;" << endl;
		 fout << "DIMENSIONS ntax=" << m_iNumber << ";" << endl;
		 fout << "TAXLABELS" << endl;
		 
		 unsigned int ix=1;
		 vector<string>::iterator iter1;
		 for(iter1 = m_gName.begin(); iter1 != m_gName.end(); ++iter1)
         {
		     fout << "[" << ix << "] " << "'" << *iter1 << "'" << endl;
			 ++ix;
		 }
		 fout << ";" << endl;
		 fout << "END; [Taxa]" << endl;
		 
		 fout << endl;
		 
		 
		 // retrive threshold
		 double dSum = 0;
         unsigned int iSplitSize = m_gSplit.size();
         
		 for(unsigned int i = iSplitSize-1; i >= iSplitSize-m_iNumber; --i)
		      dSum += m_gSplit[i].GetWeight();
		      
         double dThresh = dSum / m_iNumber;
         dThresh = dThresh * dPercent / 100;
         
         cout << dPercent <<"% of the average weight:"<< dThresh << endl;
		 
		 
		 // Output CSplits and their weights
		 fout << "BEGIN splits;" << endl;
		 fout << "DIMENSIONS ntax=" << m_iNumber << " nsplits=" << m_gSplit.size() << ";" << endl;
		 fout << "FORMAT labels=no weights=yes confidences=no intervals=no;" << endl;
		 fout << "MATRIX" << endl;
		 
		 unsigned int ia=1;
		 
		 // nontrivail splits
         for(vector<CSplit>::iterator iter = m_gSplit.begin(); iter != m_gSplit.end()-m_iNumber; ++iter)
         {
            double tWeight = (*iter).GetWeight(); 
            
            if ( tWeight >= dThresh)
            { 
                 vector<unsigned int> gLeft, gRight;
                 gLeft = (*iter).GetLblk();
	       	     gRight = (*iter).GetRblk();
		         unsigned int iLlen = gLeft.size();
			     unsigned int iRlen = gRight.size();
		         unsigned int iLen = min(iLlen, iRlen);
			
                 fout << "[" << ia << ", " << "size=" << iLen << "]" << "\t" << tWeight << "\t";
  	             fout <<"\t";
			
			
		       	if(gLeft[gLeft.size()-1] == 1)
		    	{
			        for(vector<unsigned int>::iterator it = gLeft.end()-1; it != gLeft.begin(); --it)
                        fout << *it << " ";
		             fout << gLeft[0] << "," << endl;
			     }
			    else
			    {
                    for(vector<unsigned int>::iterator it = gRight.end()-1; it != gRight.begin(); --it)
                       fout << *it << " ";
			        fout << gRight[0] << "," << endl;			
			     }
			
			    ++ia;
           }
         }
         
         // trvial splits
         for(vector<CSplit>::iterator iter = m_gSplit.end()-m_iNumber; iter != m_gSplit.end(); ++iter)
         {
            
		    vector<unsigned int> gLeft, gRight;
            gLeft = (*iter).GetLblk();
			gRight = (*iter).GetRblk();
		    unsigned int iLlen = gLeft.size();
			unsigned int iRlen = gRight.size();
			unsigned int iLen = min(iLlen, iRlen);
			
			fout << "[" << ia << ", " << "size=" << iLen << "]" << "\t" << (*iter).GetWeight() << "\t";
			fout <<"\t";
			
			
			if(gLeft[gLeft.size()-1] == 1)
			{
			     for(vector<unsigned int>::iterator it = gLeft.end()-1; it != gLeft.begin(); --it)
                     fout << *it << " ";
				  fout << gLeft[0] << "," << endl;
			}
			else
			{
              for(vector<unsigned int>::iterator it = gRight.end()-1; it != gRight.begin(); --it)
                     fout << *it << " ";
			  fout << gRight[0] << "," << endl;			
			}
			
			++ia;
         }
		   		
		 fout << ";" << endl;			
		 fout << "END; [splits]" << endl;
     }
     
}

// generate complete splits and output to file by Minimum method from Multiple alignment file
void CSplitSys::ProcessAlnMin(string& InputFile, string& OutputFile, double dPercent)
{
   // Read data from the weight file InputFile
	
   cout << "-----------------------------------------------" << endl;
   cout << "Reading Multiple alignment file starts!" << endl; 
   
   ReadAlignment(InputFile); 
   
   cout << "Reading Multiple alignment file finished!" << endl; 
   cout << "-----------------------------------------------" << endl;
   cout << endl;

   cout << "-----------------------------------------------" << endl;
   cout << "Generating quintets starts!" << endl;
   
   GenerateQuintet();
      
   cout << "Generating quintets ended!" << endl;
   cout << "-----------------------------------------------" << endl;
   cout << endl;
   
   cout << "-----------------------------------------------" << endl;
   cout << "Generating septets starts!" << endl;
   
   GenerateSeptet();
      
   cout << "Generating septets ended!" << endl;
   cout << "-----------------------------------------------" << endl;
   cout << endl;
   
   cout << "-----------------------------------------------" << endl;
   cout << "Calculating complete non-trival splits weights starts!" << endl;  
   
   ExpandAlignmentMin();
   
   cout << "Calculating complete non-trivial splits ended!" << endl;
   cout << "-----------------------------------------------" << endl;
   cout << endl;
   
   cout << "-----------------------------------------------" << endl;
   cout << "Generating trivial splits starts!" << endl; 
   
   GenerateTrivialSplit();
   
   cout << "Generating trivial splits ended!" << endl;	 
   cout << "------------------------------------------------------" << endl;
	
   cout << "------------------------------------------------------" << endl;
   cout << "Writing the complete splits to OutputFile starts!" << endl;	
	
   Output(OutputFile, dPercent);

   cout << "Writing the complete splits to OutputFile ended!" << endl;	 
   cout << "------------------------------------------------------" << endl;
}                                      



