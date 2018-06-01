#ifndef _CSPLIT_H
#define _CSPLIT_H
#include <vector>

using namespace std;

// Declaration of the class SPLIT

class CSplit{
public:

    // Default constructor
	CSplit():m_gLblk(),m_gRblk(), m_dWeight(){};       
	
	// Constructor with parameters
	CSplit(const vector<unsigned int>& gLb, const vector<unsigned int>& gRb, double dWei): m_gLblk(gLb), m_gRblk(gRb), m_dWeight(dWei){};                                   	
	
    // Copy constructor
	CSplit(const CSplit& rhs): m_gLblk(rhs.m_gLblk), m_gRblk(rhs.m_gRblk), m_dWeight(rhs.m_dWeight){};                                 
		   
    // Operator constructor
	CSplit& operator=(const CSplit& rhs)
	{
	    m_gLblk = rhs.m_gLblk;
	    m_gRblk = rhs.m_gRblk;
	    m_dWeight=rhs.m_dWeight;
	    return *this;
    }
                                                         
    vector<unsigned int> GetLblk() { return m_gLblk; }                                 // Return the left block 
    vector<unsigned int> GetRblk() { return m_gRblk; }                                 // Return the right block  
    double GetWeight()  { return m_dWeight; }                                          // Return split_weight
	
	void SetLabel(const vector<unsigned int>& gLb, const vector<unsigned int>& gRb)    // Set split label 
    { 
         m_gLblk = gLb; 
         m_gRblk = gRb;
    }   
    
    void SetWeight(double dWei){ m_dWeight = dWei;}                                    // Set split weights to be dWei
	
private:
	vector<unsigned int> m_gLblk;                                                      // The left block of the label
	vector<unsigned int> m_gRblk;                                                      // The right block of the label
	double m_dWeight;                                                                  // Split weight
};

#endif
//EOF
