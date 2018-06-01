
# QS-Net
QS-Net is a more accurate method for reconstructing phylogenetic networks extending our previous method Quartet based on Sextet. 
QS-Net requires that C++ is installed and configured on the system.

QS-Net takes Multiple Sequence Alignments(MSAs) as an input file and the .NEXU file format is used to store all splits. We can use SplitsTree to visualize the reconstructed result saved in .NEXU file.

Usage: command-line, run the following command:
       QS-Net [parameter] inputfile [parameter] threshold [parameter] outputfile
       
  example: QS-Net -f data/five1.fas -c 1 -o output/five1.fas
  
  parameters:
       -f: InputFile Name: followed by the name of the file.
       -c: Threshold to determine the splits to be filtered.
			 -o: OutputFile Name-the reconstructed result saved in it.
       

  
