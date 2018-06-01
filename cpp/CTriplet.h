#ifndef _CTRIPLET_H
#define _CTRIPLET_H

// Declaration of the class CTRIPLET

#include<vector>

using namespace std;

class CTriplet
{
public:
	CTriplet(){};                                                    // Default constructor
	CTriplet(unsigned int, unsigned int[2], double);                 // Constructor with parameters
	CTriplet(const CTriplet& );                                      // Copy constructor
	CTriplet& operator=(const CTriplet& );                           // Operator constructor
       
    const unsigned int GetLblk() const { return m_gLblk; }           // Return left block of the triplet label
	const unsigned int* GetRblk() const { return m_gRblk; }          // Return right block of the triplet label
	double GetWeight()  { return m_dWeight; }                        // Return triplet weight
	void SetLabel(unsigned int, unsigned int[2]);                    // Set triplet label 
    void SetWeight(double dWei){ m_dWeight = dWei;}                  // Set split weight to be dWei
	void DecreaseOrder();                                            // Ordering according to the decreasing order 2|14 to 2|41
	unsigned int GetIndex() const;                                   // Return the index
       
private:   
	unsigned int m_gLblk;                                             // The left block of the label
	unsigned int m_gRblk[2];                                          // The right block of the label
	double m_dWeight;                                                 // Triplet weight
};

#endif
// EOF
