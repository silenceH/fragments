CFLAGS=-Wall -g

frags:
	g++ -std=c++11 -o frags_11 frags_11.cpp -I$(RDBASE)/Code -I$(RDBASE)/EXtern -L$(RDBASE)/lib -lFileParsers -lGraphMol -lRDGeneral -lChemTransforms
# make the test file to compile the Fragment analysis

kenn:
	g++ -std=c++11 -o kenn kennewell_implementation.cpp -I$(RDBASE)/Code -I$(RDBASE)/EXtern -L$(RDBASE)/lib -lFileParsers -lGraphMol -lRDGeneral -lChemTransforms
# need to make an object file Fragments.o 
Fragments.o : Fragments.cpp Fragments.hpp
	g++ -std=c++11 -c -O Fragments.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lFileParsers -lGraphMol -lRDGeneral -lChemTransforms

# Then we need to compile the test code
Test_Fragments.o : Test_Fragments.cpp Fragments.hpp 
	g++ -std=c++11 -c -O Test_Fragments.cpp

Test_Fragments : Fragments.o Test_Fragments.o
	g++ -O -o Test_Fragments Fragments.o Test_Fragments.o
	
clean:
	rm -f kenn frags_11 test file_test
