# Project: QS-Net
# Makefile created by Dev-C++ 5.9.2

CPP      = /data2/tools/gcc-6.2.0/makehere/bin/g++
CC       = /data2/tools/gcc-6.2.0/makehere/bin/gcc
OBJ      = cpp/auxilary.o cpp/CalWeiMin.o cpp/CQuartet.o cpp/CQuintet.o cpp/CSeptet.o cpp/CSextet.o cpp/CSplitSys.o cpp/CTriplet.o cpp/ExpandAlignmentMin.o cpp/GenerateQuintet.o cpp/GenerateSeptet.o cpp/GenerateTrivialSplit.o cpp/main.o cpp/MinQuintet.o cpp/MinSeptet.o cpp/ReadAlignment.o
LINKOBJ  = cpp/auxilary.o cpp/CalWeiMin.o cpp/CQuartet.o cpp/CQuintet.o cpp/CSeptet.o cpp/CSextet.o cpp/CSplitSys.o cpp/CTriplet.o cpp/ExpandAlignmentMin.o cpp/GenerateQuintet.o cpp/GenerateSeptet.o cpp/GenerateTrivialSplit.o cpp/main.o cpp/MinQuintet.o cpp/MinSeptet.o cpp/ReadAlignment.o
LIBS     = -L"/data2/tools/" -static-libgcc
INCS     = -I"/data2/tools/"
CXXINCS  = -I"/data2/tools/"
BIN      = QS-Net
CXXFLAGS = $(CXXINCS) 
CFLAGS   = $(INCS) 
RM       = rm -f

.PHONY: all all-before all-after clean clean-custom

all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(LINKOBJ) -o $(BIN) $(LIBS)

cpp/auxilary.o: cpp/auxilary.cpp
	$(CPP) -c cpp/auxilary.cpp -o cpp/auxilary.o $(CXXFLAGS)

cpp/CalWeiMin.o: cpp/CalWeiMin.cpp
	$(CPP) -c cpp/CalWeiMin.cpp -o cpp/CalWeiMin.o $(CXXFLAGS)

cpp/CQuartet.o: cpp/CQuartet.cpp
	$(CPP) -c cpp/CQuartet.cpp -o cpp/CQuartet.o $(CXXFLAGS)

cpp/CQuintet.o: cpp/CQuintet.cpp
	$(CPP) -c cpp/CQuintet.cpp -o cpp/CQuintet.o $(CXXFLAGS)

cpp/CSeptet.o: cpp/CSeptet.cpp
	$(CPP) -c cpp/CSeptet.cpp -o cpp/CSeptet.o $(CXXFLAGS)

cpp/CSextet.o: cpp/CSextet.cpp
	$(CPP) -c cpp/CSextet.cpp -o cpp/CSextet.o $(CXXFLAGS)

cpp/CSplitSys.o: cpp/CSplitSys.cpp
	$(CPP) -c cpp/CSplitSys.cpp -o cpp/CSplitSys.o $(CXXFLAGS)

cpp/CTriplet.o: cpp/CTriplet.cpp
	$(CPP) -c cpp/CTriplet.cpp -o cpp/CTriplet.o $(CXXFLAGS)

cpp/ExpandAlignmentMin.o: cpp/ExpandAlignmentMin.cpp
	$(CPP) -c cpp/ExpandAlignmentMin.cpp -o cpp/ExpandAlignmentMin.o $(CXXFLAGS)

cpp/GenerateQuintet.o: cpp/GenerateQuintet.cpp
	$(CPP) -c cpp/GenerateQuintet.cpp -o cpp/GenerateQuintet.o $(CXXFLAGS)

cpp/GenerateSeptet.o: cpp/GenerateSeptet.cpp
	$(CPP) -c cpp/GenerateSeptet.cpp -o cpp/GenerateSeptet.o $(CXXFLAGS)

cpp/GenerateTrivialSplit.o: cpp/GenerateTrivialSplit.cpp
	$(CPP) -c cpp/GenerateTrivialSplit.cpp -o cpp/GenerateTrivialSplit.o $(CXXFLAGS)

cpp/main.o: cpp/main.cpp
	$(CPP) -c cpp/main.cpp -o cpp/main.o $(CXXFLAGS)

cpp/MinQuintet.o: cpp/MinQuintet.cpp
	$(CPP) -c cpp/MinQuintet.cpp -o cpp/MinQuintet.o $(CXXFLAGS)

cpp/MinSeptet.o: cpp/MinSeptet.cpp
	$(CPP) -c cpp/MinSeptet.cpp -o cpp/MinSeptet.o $(CXXFLAGS)

cpp/ReadAlignment.o: cpp/ReadAlignment.cpp
	$(CPP) -c cpp/ReadAlignment.cpp -o cpp/ReadAlignment.o $(CXXFLAGS)
