#ifndef _CQUARTET_H
#define _CQUARTET_H

// Declaration of the class QUARTET

class CQuartet
{
public:
	CQuartet(){}                                                     // Default constructor
	CQuartet(unsigned int[2], unsigned int[2], double);              // Constructor with parameters
	CQuartet(const CQuartet& );                                      // Copy constructor
	CQuartet& operator=(const CQuartet& );                           // Operator constructor
       
    const unsigned int* GetLblk() const { return m_gLblk; }          // Return left block of the quartet label
	const unsigned int* GetRblk() const { return m_gRblk; }          // Return right block of the quartet label
	double GetWeight()  { return m_dWeight; }                        // Return quartet_weight
	void SetLabel(unsigned int[2], unsigned int[2]);                 // Set quartet label 
    void SetWeight(double dWei){ m_dWeight = dWei;}                  // Set split weights to be dWei
	void DecreaseOrder();                                            // Ordering according to the decreasing order 32|14 to 41|32
	unsigned int GetIndex() const;                                   // Return the index
       
private:   
	unsigned int m_gLblk[2];                                         // The left block of the label
	unsigned int m_gRblk[2];                                         // The right block of the label
	double m_dWeight;                                                // Quartet weight
};

#endif
// EOF
