#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib> 
#include <vector>

#include "CSplitSys.h"
#include "auxilary.h"

using namespace std;

void CSplitSys::ReadWeight(const string& sInputFile)
{   
      cout << "----------------------------------------"<< endl;
      cout << "Reading Weight File Starts!" << endl;  
      
     // Open file sInputFile for reading
     ifstream fin(sInputFile.c_str());                                 
               
     if(!fin ) 
     {   
         cerr << "Error to  open " << sInputFile.c_str() << " for reading!" << endl;  
         system("pause"); 
         exit(-1);  
     }
     else
     {  
        string sLine, sWord, sTemp;           // sLine stores a line. sWord stores words. 
       
        // Define the identifiers
        const string cTAXANUMBER = "taxanumber:"; 
        const string cTAXON = "taxon:";
        const string cQUARTET = "quartet:";   
        const string cTRIPLET = "triplet:";
		       
         // Read the file line by line        
         while( getline(fin,sLine) )                                         
         {   
               if(sLine == "") continue;    // if it is a blank line, go back and read the next line
			   
               // Read the line words by words and store them to a string vector words.
               istringstream stream(sLine);                                   
               vector<string> gWords;
               unsigned int iSize;
                    
               while(stream >> sWord)
               {
                    gWords.push_back(sWord.c_str());
               }
                                          
              // Case 1: the line cotains informaiton of splits, read split name and split weight. 
			  //         e.g. quartet: 001 002 003 004 weights: 0.1756 1.2800223125448445E-71 3.495951593252628E-77;
			  
              if(gWords[0] == cQUARTET)                              
               {
                      unsigned int gLeft[2], gRight[2];   // stores three types of quartets 41|32, 42|31, 43|21.
                      
                     // Break into 3 types of quartets and their weight 4321: 41|32, 42|31, 43|21.           
                     // The first quarte 43|21 and its weight
                     gLeft[0] = atoi(gWords[1].c_str());
                     gLeft[1] = atoi(gWords[2].c_str());
                     gRight[0] = atoi(gWords[3].c_str());
                     gRight[1] = atoi(gWords[4].c_str());
                     
                     CQuartet qOne(gLeft, gRight, atof(gWords[6].c_str()));
                     m_gQuartet.push_back(qOne);
                     
                     // The second quarte 42|31 and its weight
                     swap(gLeft[1], gRight[0]);
                     CQuartet qTwo(gLeft, gRight, atof(gWords[7].c_str()));
                     m_gQuartet.push_back(qTwo);
                                  
                     // The third quarte 41|32 and its weight
                     swap(gLeft[1], gRight[1]);
                     sTemp = gWords[8];
                     iSize = sTemp.size()-1; 
                     if( sTemp[iSize] == ';' )
					     sTemp = sTemp.substr(0, iSize);
                     CQuartet qThree(gLeft, gRight, atof(sTemp.c_str()));
                     m_gQuartet.push_back(qThree);                                                      
              }   
              
			  // Case 2: the line cotains informaiton of splits, read split name and split weight. 
			  //         e.g. triplet: 001 002 003  weights: 0.1756 1.2800223125448445E-71 3.495951593252628E-77;
			  
              else if(gWords[0] == cTRIPLET)                              
               {
                      unsigned int gLeft, gRight[2];   // stores three types of triplets 1|32, 2|31, 3|21.
                      
                     // Break into 3 types of triplets and their weight 321: 1|32, 2|31, 3|21.           
                     // The first triplet 3|21 and its weight
                     gLeft = atoi(gWords[1].c_str());
                     gRight[0] = atoi(gWords[2].c_str());
                     gRight[1] = atoi(gWords[3].c_str());
                     
                     CTriplet tOne(gLeft, gRight, atof(gWords[5].c_str()));
                     m_gTriplet.push_back(tOne);
                     
                     // The second triplet 2|31 and its weight
                     swap(gLeft, gRight[0]);
                     CTriplet tTwo(gLeft, gRight, atof(gWords[6].c_str()));
                     m_gTriplet.push_back(tTwo);
                                  
                     // The third quartes 1|32 and its weight
                     swap(gLeft, gRight[1]);
                     sTemp = gWords[7];
                     iSize = sTemp.size()-1; 
                     if( sTemp[iSize] == ';' )
					     sTemp = sTemp.substr(0, iSize);
                     CTriplet tThree(gLeft, gRight, atof(sTemp.c_str()));
                     m_gTriplet.push_back(tThree);                                                      
              }   

              // Case 3: the line cotains informaiton of taxon and its name. 
			  //       e.g. taxon:   003   name: NC_005958 ;
              else if( gWords[0] == cTAXON )
			  {
                   string sName= gWords[3];                          // name of the taxon
                   iSize = sName.size()-1; 
				   if(sName[iSize] == ';') sName = sName.substr(0, iSize);
                   m_gName.push_back(sName);                                       
			   }               
              // Case 4: the line cotains informaiton of taxa number.
			  //        e.g. taxanumber: 9;
			   else if (gWords[0] == cTAXANUMBER)
			   {
					sTemp = gWords[1].substr(0, gWords[1].size()-1);
					m_iNumber = atoi(sTemp.c_str());
				}           
         }
     }
     fin.close();   
     
     // Reorder quartets 
     OrderQuartets();
	 
	 if(m_gTriplet.size() == 0)
	 {
	     string::size_type i, j, k; 
		 unsigned int gRight[2];
		 double dWei = 0;
		 
		for(i = 1; i != m_iNumber-1; ++i)
		   for( j = i+1; j != m_iNumber; ++j)
			   for( k = j+1; k != m_iNumber +1; ++k)
		{
		      // initialize 3 triplets 
			  
			  // Triplet 1: i|kj, e.g. 1|32
			  gRight[0] = k;
			  gRight[1] = j;
			  
		      CTriplet tOne(i, gRight, dWei);
		      
		      // push triplet 1 to the corresponding position of  m_gTriplet
			  m_gTriplet.push_back(tOne);
			  
			  // push triplet 2: j|ki to m_gTriplet 2|31
              gRight[1] = i;
              
		      CTriplet tTwo(j, gRight, dWei);
			  m_gTriplet.push_back(tTwo);
			  
			  // push triplet 3: k|ji to m_gTriplet 3|21
              gRight[0] = j;
              
		      CTriplet tThree(k, gRight, dWei);
			  m_gTriplet.push_back(tThree);
		}

	 }
	
	 // reorder triplets 
	 OrderTriplets();
     
     cout<< "Reading Weight File End!" << endl;         
     cout << "----------------------------------------"<< endl;
     
}
          
//EOF
	  
	 
     

