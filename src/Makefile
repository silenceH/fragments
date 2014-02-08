CFLAGS=-Wall -g

# make the test file to compile the Fragment analysis
# need to make an object file Fragments.o 
Fragments.o : Fragments.cpp Fragments.hpp
	g++ -std=c++11 -c -O Fragments.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lFileParsers -lGraphMol -lRDGeneral -lChemTransforms

# Then we need to compile the test code
Test_Fragments.o : Test_Fragments.cpp Fragments.hpp 
	g++ -std=c++11 -c -O Test_Fragments.cpp

Test_Fragments : Fragments.o Test_Fragments.o
	g++ -O -o Test_Fragments Fragments.o Test_Fragments.o
	
