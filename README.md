
# QS-Net

User Manual for QS-Net
Introduction
The QS-Net is a phylogenetic network reconstruction method taking advantage of information on the relationship among six taxa. 

Compile and run 
The program runs under Windows. The user can directly run the executable file "QS-Net.exe" in the command window. All source code is in the "cpp" folder, which has been compiled and can be run directly. And the software used is Dev-cpp, which can be downloaded and installed on the website https://sourceforge.net/projects/orwelldevcpp/ , or you can search for other versions directly by browser.

Use from command-line
QS-Net takes Multiple Sequence Alignments (MSAs) as an input file and the .NEXU file format is used to store all splits.
Command description:
QS-Net  [-PARAMS PARAM_VALUE] 
PARAMS:
        -f   InputFile Name
        -c   Threshold to determine the splits to be filtered.
        -o   Output File Name 
Example:  QS-Net -f data/five1.fas -c 1 -o output/five1.fas

Running steps
First download the code and data on the github website https://github.com/Tmyiri/QS-Net and save it on a disk (for example, directly on the desktop). Next open a command line window, enter the command to jump to the specified area of the workspace, as shown below:

Next, enter the command in the example.

After the program is executed, there will be a five1.nex file in the output folder, and the result of the build will be saved in this file.
Visualization
The output file will be a file in Nexus format, which can be viewed by SplitsTree http://www.splitstree.org/. 
  
