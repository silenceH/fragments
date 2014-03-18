CFLAGS=-Wall -g

#frags:
#	g++ -std=c++11 -o frags_11 frags_11.cpp -I$(RDBASE)/Code -I$(RDBASE)/EXtern -L$(RDBASE)/lib -lFileParsers -lGraphMol -lRDGeneral -lChemTransforms
## make the test file to compile the Fragment analysis
#
#kenn:
#	g++ -std=c++11 -o kenn kennewell_implementation.cpp -I$(RDBASE)/Code -I$(RDBASE)/EXtern -L$(RDBASE)/lib -lFileParsers -lGraphMol -lRDGeneral -lChemTransforms
# need to make an object file Fragments.o 
Frag.o : Frag.cpp Frag.hpp
	g++ -std=c++11 -c -O Frag.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lChemTransforms

RDKMols.o : RDKMols.cpp RDKMols.hpp
	g++ -std=c++11 -c -O RDKMols.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lFileParsers -lChemTransforms

Bioisosteres.o: Bioisosteres.cpp Bioisosteres.hpp
	g++ -std=c++11 -c -O Bioisosteres.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lFileParsers -lChemTransforms

FragGroup.o: FragGroup.cpp FragGroup.hpp
	g++ -std=c++11 -c -O FragGroup.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lFileParsers -lChemTransforms

# Then we need to compile the test code
main.o : main.cpp Frag.hpp RDKMols.hpp Bioisosteres.hpp FragGroup.hpp
	g++ -std=c++11 -c -O main.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lFileParsers -lChemTransforms

main: Frag.o RDKMols.o Bioisosteres.o FragGroup.o main.o
	g++ -O -o main Frag.o Bioisosteres.o RDKMols.o FragGroup.o main.o -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lChemTransforms -lFileParsers 
	
clean:
	rm -f kenn frags_11 test file_test main
