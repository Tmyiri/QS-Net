#ifndef _CQUINTET_H
#define _CQUINTET_H

// Declaration of the class CQuintet

class CQuintet
{
public:
	CQuintet(){}                                                     // Default constructor
	CQuintet(unsigned int[2], unsigned int[3], double);              // Constructor with parameters
	CQuintet(const CQuintet& );                                      // Copy constructor
	CQuintet& operator=(const CQuintet& );                           // Operator constructor
       
    unsigned int* GetLblk() { return m_gLblk; }                      // Return left block of the quintet label
	unsigned int* GetRblk() { return m_gRblk; }                      // Return right block of the quintet label
	double GetWeight()  { return m_dWeight; }                        // Return quintet_weight
	void SetLabel(unsigned int[2], unsigned int[3]);                 // Set split label 
    void SetWeight(double dWei){ m_dWeight = dWei;}                  // Set split weights to be dWei
	void DecreaseOrder();                                            // Ordering according to the decreasing order 23|514 to 32|541
	unsigned int GetIndex() const;                                   // Return the index
       
private:   
	unsigned int m_gLblk[2];                                         // The left block of the label
	unsigned int m_gRblk[3];                                         // The right block of the label
	double m_dWeight;                                                // Split weight
};

#endif
// EOF
