CFLAGS=-Wall -g

frags:
	g++ -std=c++11 -o frags_11 frags_11.cpp -I$(RDBASE)/Code -I$(RDBASE)/EXtern -L$(RDBASE)/lib \
		-lFileParsers -lGraphMol -lRDGeneral -lChemTransforms -lFingerprints -lDataStructs
## make the test file to compile the Fragment analysis
#
#kenn:
#	g++ -std=c++11 -o kenn kennewell_implementation.cpp -I$(RDBASE)/Code -I$(RDBASE)/EXtern -L$(RDBASE)/lib -lFileParsers -lGraphMol -lRDGeneral -lChemTransforms

bioisosteres :
	g++ -std=c++11 -o bioisosteres bioisosteres.cpp -I$(RDBASE)/Code -I$(RDBASE)/EXtern -L$(RDBASE)/lib \
		-lFileParsers -lGraphMol -lRDGeneral -lChemTransforms -lFingerprints -lDataStructs
	
# need to make an object file Fragments.o 
Frag.o : Frag.cpp Frag.hpp
	g++ -std=c++11 -c -O Frag.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lChemTransforms -lFingerprints

Bioisosteres.o: Bioisosteres.cpp Bioisosteres.hpp
	g++ -std=c++11 -c -O Bioisosteres.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lFileParsers -lChemTransforms

FragGroup.o: FragGroup.cpp FragGroup.hpp
	g++ -std=c++11 -c -O FragGroup.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lFileParsers -lChemTransforms

# Then we need to compile the test code
main.o : main.cpp Frag.hpp Bioisosteres.hpp FragGroup.hpp
	g++ -std=c++11 -c -O main.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lFileParsers -lChemTransforms -lFingerprints

main: Frag.o Bioisosteres.o FragGroup.o main.o
	g++ -O -o main Frag.o Bioisosteres.o FragGroup.o main.o -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lChemTransforms -lFileParsers -lFingerprints
	
tanimoto.o: tanimoto.cpp Frag.hpp Bioisosteres.hpp
	g++ -std=c++11 -c -O tanimoto.cpp -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lFileParsers -lChemTransforms -lFingerprints

tanimoto: Frag.o Bioisosteres.o FragGroup.o main.o
	g++ -O -o tanimoto Frag.o Bioisosteres.o tanimoto.o -I$(RDBASE)/Code -I$(RDBASE)/Extern -L$(RDBASE)/lib -lGraphMol -lRDGeneral -lChemTransforms -lFileParsers -lFingerprints

clean:
	rm -f kenn frags_11 test file_test main tanimoto bioisosteres
