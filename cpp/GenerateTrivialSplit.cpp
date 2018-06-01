#include <string>
#include <iostream>
#include "CSplitSys.h"
#include "CTriplet.h"

using namespace std;

// Gerneate all trivial splits from triplets and non-trivial full splits
void CSplitSys::GenerateTrivialSplit()
{
	CTriplet tEmpty;
	unsigned int gTemp[2];
		
	vector<double> dWeight;
	
	// generate the weight vector of m_gTriplet by w(a | X-a) = min_{b,c \in X-a} {w(a|bc)-\sum_{a \in A, bc \in B}A|B}
	for(vector<CTriplet>::iterator iter = m_gTriplet.begin(); iter != m_gTriplet.end(); ++iter)
	     dWeight.push_back((*iter).GetWeight());
	
	// go throught all the full splits and minus the corresponding weights from w(a|bc)
	for(vector<CSplit>::iterator iter = m_gSplit.begin(); iter != m_gSplit.end(); ++iter)
	{
	    vector<unsigned int> gLeft = (*iter).GetLblk(); 
		vector<unsigned int> gRight = (*iter).GetRblk();
		double dWei = (*iter).GetWeight();
				
	    for(vector<unsigned int>::iterator iter1= gLeft.begin(); iter1!= gLeft.end(); ++iter1)
		    for(vector<unsigned int>::iterator iter2 = gRight.begin(); iter2 != gRight.end()-1; ++iter2)
			   for(vector<unsigned int>::iterator iter3= iter2+1; iter3 != gRight.end(); ++iter3)
		{						
		    gTemp[0] = *(iter2);
			gTemp[1] = *(iter3);
		    tEmpty.SetLabel(*(iter1), gTemp);
			dWeight[tEmpty.GetIndex()] -= dWei;
		}
		
		
		for(vector<unsigned int>::iterator iter1= gRight.begin(); iter1!= gRight.end(); ++iter1)
		    for(vector<unsigned int>::iterator iter2= gLeft.begin(); iter2!= gLeft.end()-1; ++iter2)
			   for(vector<unsigned int>::iterator iter3= iter2+1; iter3!= gLeft.end(); ++iter3)
		{
			gTemp[0] = *(iter2);
			gTemp[1] = *(iter3);
		    tEmpty.SetLabel(*(iter1), gTemp);
			
			dWeight[tEmpty.GetIndex()] -= dWei;
		}
	}
	
	// generating all the trivial full splits
	for(unsigned int i = 1; i != m_iNumber+1; ++i)
	{
	    double dWei = 10000000;
		vector<unsigned int> gLeft, gRight;
                
		gLeft.push_back(i);
	            
		for(unsigned j = 1; j != m_iNumber+1; ++j)
		{
		    if (j == i) continue;
			else
			{
				for(unsigned k = j+1; k != m_iNumber+1; ++k)
				{
					if( k == i) continue;
					else
					  {
					      gTemp[0] = j;
						  gTemp[1] = k;
						  
						  tEmpty.SetLabel(i, gTemp);
						  
						  double dTemp =  dWeight[tEmpty.GetIndex()];
						  
						  if (dTemp <= 0) 
						  {
						      dWei =1;
							  break;
						  }
						  else
						  {
						      if(dWei > dTemp) dWei = dTemp;
						  }
					  }
				}
			}
                        
			gRight.push_back(j);
		}
                  
		// w(1|234...n) = min_{j,k \in 2,3..n} w(1|jk)
		CSplit sTemp(gLeft, gRight, dWei);
               m_gSplit.push_back(sTemp);
        }
}
