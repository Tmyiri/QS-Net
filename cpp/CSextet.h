#ifndef _CSEXTET_H
#define _CSEXTET_H

#include<vector>

using namespace std;

class CSextet
{
public:
	CSextet(){};                                                    // Default constructor
	CSextet(unsigned int[3], unsigned int[3], double);              // Constructor with parameters
	CSextet(const CSextet& );                                       // Copy constructor
	CSextet& operator=(const CSextet& );                            // Operator constructor
       
    const unsigned int* GetLblk() const { return m_gLblk; }         // Return left block of the sextet label
	const unsigned int* GetRblk() const { return m_gRblk; }         // Return right block of the sextet label
	double GetWeight()  { return m_dWeight; }                       // Return sextet weight
	void SetLabel(unsigned int[3], unsigned int[3]);                // Set sextet label 
    void SetWeight(double dWei){ m_dWeight = dWei;}                 // Set split weights to be dWei
	void DecreaseOrder();                                           // Ordering according to the decreasing order 123|456 to 654|321
	unsigned int GetIndex() const;                                  // Return the index
       
private:   
	unsigned int m_gLblk[3];                                        // The left block of the label
	unsigned int m_gRblk[3];                                        // The right block of the label
	double m_dWeight;                                               // sextet weight
};

#endif
