The method of running QS-Net on a Linux system is the same as the Windows system. 

We extended the Makefile and used the "make" command to generate a  "QS-Net" file. 

Then，the user can run QS-Net as follows:
   Command description: 
           QS-Net [-PARAMS PARAM_VALUE] PARAMS: 
                          -f InputFile Name 
                          -c Threshold to determine the splits to be filtered
                          -o Output File Name 
           Example: QS-Net -f data/five1.fas -c 1 -o output/five1.fas
--------------------------------------------------------------------------------------------
* If you want to compile the makefile yourself on more systems (such as Mac), you need to install G++ to compile C++ code. Then the next step is to change the following parameters in the Makefile , finally, use the "make" command to generate a  "QS-Net" file.
       CPP      = /data2/tools/gcc-6.2.0/makehere/bin/g++


       CC       = /data2/tools/gcc-6.2.0/makehere/bin/gcc


       LIBS     = -L"/data2/tools/" -static-libgcc
       

INCS     = -I"/data2/tools/"

       CXXINCS  = -I"/data2/tools/"