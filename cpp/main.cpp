#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "CSplitSys.h"

using namespace std;

int main(unsigned int argc, char *argv[])
{
    string InputFile = "";
    string OutputFile = "output";
    string pStr;
    double Threshold = 0.01;
    
    cout<< argc << endl;
    
    if(argc % 2 == 0)
    {       cout << "Usage: QuartetMethods [-params values] " << endl;
            cout << "         -f File Name: followed by the name of the file" << endl;
            cout << "         -c Threshold to determine the splits to be filtered" << endl;
			cout << "         -o Output File Name" << endl;
            exit(0);
    }
    else
    {
        // Read options from command line:
            for(unsigned int n = 1; n < argc; n = n + 2)
            {
                if(argv[n][0] != '-')
                {
                   cout << " Parameter must start with -" << endl;
                   exit(0);
                }
                else
                {	
                   switch(argv[n][1])
                   {
                       case 'f': // Read file name
                            InputFile = argv[n+1];
		                    break;
                       case 'o':  // output file name
                            OutputFile = argv[n+1];
                            break; 		
                       case 'c':  // Read Threshold
                            pStr = argv[n+1];
                            Threshold = atof(pStr.c_str());
                            break;	
   			        	default: // No such option
					        cout << "Illegal option: " << argv[n] << " ! " << endl;
					        exit(0);
                   }
                }
		    }
    }
	
	
    // Codes start here
	
    // Record running time
    clock_t cStart, cFinish;
    double dTotalTime;
    cStart = clock();
	
	CSplitSys QS_Sys;
	
	if (InputFile == "")
    {   
        cout << "Please provide a file name with parameter -f!" << endl;
        exit(0);
    }
        
    QS_Sys.ProcessAlnMin(InputFile, OutputFile, Threshold);      
     
    cFinish  = clock();
    
    dTotalTime = (double)(cFinish-cStart)/CLOCKS_PER_SEC;
    
    cout << "The total Time is " << dTotalTime << " Seconds." << endl;
    
    return 0;
}
