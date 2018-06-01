#ifndef _CSEPTET_H
#define _CSEPTET_H

#include<vector>

using namespace std;

class CSeptet
{
public:
	CSeptet(){};                                                    // Default constructor
	CSeptet(unsigned int[3], unsigned int[4], double);              // Constructor with parameters
	CSeptet(const CSeptet& );                                       // Copy constructor
	CSeptet& operator=(const CSeptet& );                            // Operator constructor
       
    unsigned int* GetLblk()  { return m_gLblk; }                    // Return left block of the septet label
	unsigned int* GetRblk()  { return m_gRblk; }                    // Return right block of the septet label
	double GetWeight() const { return m_dWeight; }                  // Return septet weight
	void SetLabel(unsigned int[3], unsigned int[4]);                // Set septet label 
    void SetWeight(double dWei){ m_dWeight = dWei;}                 // Set split weight to be dWei
	void DecreaseOrder();                                           // Ordering according to the decreasing order 123|4567 to 321|7654
	unsigned int GetIndex() const;                                  // Return the index
       
private:   
	unsigned int m_gLblk[3];                                        // The left block of the label
	unsigned int m_gRblk[4];                                        // The right block of the label
	double m_dWeight;                                               // septet weight
};

#endif
