#ifndef _CSPLIT_SYS_H
#define _CSPLIT_SYS_H

// Body of CSplitSys.h
#include <string>
#include "CTriplet.h"
#include "CQuartet.h"
#include "CQuintet.h"
#include "CSextet.h"
#include "CSeptet.h"
#include "CSplit.h"

using namespace std;

class CSplitSys
{
public:
    // Default constructor
	CSplitSys(): m_iNumber(), m_gName(), m_gTriplet(), m_gQuartet(), m_gQuintet(), m_gSextet(), m_gSeptet(), m_gSplit(){} 
	 
    // Constructor with parameters
    CSplitSys(unsigned int iNumber, const vector<string>& gName, const vector<CTriplet>& gTriplet, 
	const vector<CQuartet>& gQuartet, const vector<CQuintet>& gQuintet,const vector<CSextet>& gSextet,const vector<CSeptet>& gSeptet, 
	const vector<CSplit>& gSplit): m_iNumber(iNumber), m_gName(gName), m_gTriplet(gTriplet), m_gQuartet(gQuartet),m_gQuintet(gQuintet),m_gSextet(gSextet),m_gSeptet(gSeptet), m_gSplit(gSplit){}
    
    // Copy constructor                    
	CSplitSys(const CSplitSys& rhs): m_iNumber(rhs.m_iNumber), m_gName(rhs.m_gName), m_gTriplet(rhs.m_gTriplet), 
	m_gQuartet(rhs.m_gQuartet), m_gQuintet(rhs.m_gQuintet), m_gSextet(rhs.m_gSextet), m_gSeptet(rhs.m_gSeptet), m_gSplit(rhs.m_gSplit){}         
                                    
	CSplitSys& operator=(const CSplitSys& );                                     // Operator constructor

    // Set and return member variables
    unsigned int GetNumber() const { return m_iNumber;}                          // Return taxa number
    vector<string> GetName() const { return m_gName;}                            // Return taxa names
	vector<CTriplet> GetTriplet() const { return m_gTriplet;}                    // Return all Triplets
    vector<CQuartet> GetQuartet() const { return m_gQuartet;}                    // Return all quartets 
    vector<CQuintet> GetQuintet() const { return m_gQuintet;}                    // Return all quintets 
	vector<CSextet> GetSextet() const { return m_gSextet;}                       // Return all sextets
	vector<CSeptet> GetSeptet() const { return m_gSeptet;}                       // Return all septets
	vector<CSplit> GetSplit() const { return m_gSplit;}                          // Return full splits     
    void SetNumber(unsigned int iNumber){ m_iNumber = iNumber;}                  // Set taxa number to be iNumber
    void SetName(const vector<string>& gName){ m_gName = gName; }                // Set taxa names to be vector gName
	void SetTriplet(const vector<CTriplet>& gTriplet) {m_gTriplet = gTriplet;}   // Set Triplets to be vector gTriplet
	void SetQuartet(const vector<CQuartet>& gQuartet) {m_gQuartet = gQuartet;}   // Set Quartets to be vector gQuartet
	void SetSextet(const vector<CSextet>& gSextet) {m_gSextet = gSextet;}        // Set Sextets  to be vector gSextet
	void SetSeptet(const vector<CSeptet>& gSeptet) {m_gSeptet = gSeptet;}        // Set Septets  to be vector gSeptet
    
    // Ordering, shrinking and expand split system 
    void OrderQuartets();                                                        // Reoder the quatets in the split system according to their index 
	void OrderTriplets();                                                        // Reoder the triplets in the split system according to their index 
	void OrderSextets();                                                         // Reoder the sextets in the split system according to their index
	void ShrinkSplits(vector<CSplit>&);                                          // Delete splits that have weight 0
    void ExpandAlignmentMin();                                                   // Generating the non-trivial complete splits
	void GenerateTrivialSplit();                                                 // Generating trivial complete splits
      
    // Import data from files
    void ReadAlignment(const string&);                                           // save data from multiple alignment file: fasta format
   
     // Weight functions
     double MinQuintet(unsigned int gLeft[2], unsigned int gRight[3]);           // calcualte weights of split types 2|3 using minimum weight function
     double MinSeptet(unsigned int gLeft[3], unsigned int gRight[4]);            // calcualte weights of split types 3|4 using minimum weight function
     double MinWei35(unsigned int gLeft[3], unsigned int gRight[5]);             // Calculate split of type 3|5 by minimun method
     double MinWei2n(const vector<unsigned int>&, const vector<unsigned int>&);  // Calculate split of type 2|n by minimum method
     double MinWei3n(const vector<unsigned int>&, const vector<unsigned int>&);  // Calculate split of type 3|n by minimum method
     double MinWei44(unsigned int gLeft[4], unsigned int gRight[4]);             // Calculate split of type 4|4 by minimum method
     double MinWeimn(const vector<unsigned int>&, const vector<unsigned int>&);  // Calculate split of type m|n by minimum method
     void GenerateQuintet();                                                     // Generate and order all 2|3 splits using minimum method
     void GenerateSeptet();                                                      // Generate and order all 3|4 splits using minimum method
     
	 // Output
	 void Output(string&, double);                                               // Output the complete splits into OutputFile.nex
	 void ProcessAlnMin(string&, string&, double);                               // generate complete splits and output to file by the minimum method from multiple alignment file
                                                                        
private:
    unsigned int m_iNumber;                                                      // Number of taxa in the weighted system
    vector<string> m_gName;                                                      // Names of ordered taxa
	vector<CTriplet> m_gTriplet;                                                 // Ordered triplet i.e. 1|2 splits
    vector<CQuartet> m_gQuartet;                                                 // Ordered quartets 2|2 splits
    vector<CQuintet> m_gQuintet;                                                 // Ordered quintets 2|3 splits
    vector<CSextet> m_gSextet;                                                   // Ordered quintets 3|3 splits
    vector<CSeptet> m_gSeptet;                                                   // Ordered quintets 3|4 splits
	vector<CSplit> m_gSplit;                                                     // Ordered complete splits
};

#endif
// EOF


