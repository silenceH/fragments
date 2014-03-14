#include "Frag.hpp"

// constructors
Frag::Frag()
{
	file = "empty";
}

Frag::Frag(sp_fragments fragment)
{
	frag = fragment;
}


Frag::Frag(sp_fragments fragment, std::string file_name)
{
	frag = fragment;
	file = file_name;
}


// set methods
void Frag::setFile(std::string file_name)
{
	file = file_name;
}

void Frag::setSmiles(std::string s)
{
	smiles = s;
}

// get methods
std::string Frag::getSmiles() 
{
	return smiles;
}

std::string Frag::getFile() 
{
	return file;
}
