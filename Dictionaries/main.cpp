#include <fstream>
#include <iostream>

#include "Dictionary.cpp"
#include "GaborDictionary.cpp"
#include "DctDictionary.cpp"
#include "SplineDictionary.cpp"
#include "MinSplineDictionary.cpp"
#include "TrigSplineDictionary.cpp"

int main(int argc, char** argv) {
	
	if (argc != 4)	{
		return -1;
	}
	
	int atomsCount =  std::stoi(argv[1]);
	int signalSize = std::stoi(argv[2]);
	char dictType = argv[3][0];
	
	double *atoms = new double[atomsCount * signalSize];
	
	switch(dictType){
		case 'd':{
			DctDictionary generator(atomsCount, signalSize);
			generator.CreateDictionary(atoms);
			break;
		}
		case 'g':{
			GaborDictionary generator(atomsCount, signalSize);
			generator.CreateDictionary(atoms);
			break;
		}
		case 's':{
			SplineDictionary generator(atomsCount, signalSize);
			generator.CreateDictionary(atoms);
			break;
		}
		case 'm':{
			MinSplineDictionary generator(atomsCount, signalSize);
			generator.CreateDictionary(atoms);
			break;
		}
		case 't':{
			TrigSplineDictionary generator(atomsCount, signalSize);
			generator.CreateDictionary(atoms);
			break;
		}
		default:{
			std::cout << "dictType is not defined";
			return 1;
		}
	}
	
	auto fname = dictType + std::to_string(atomsCount) + std::to_string(signalSize) + ".txt";
	std::ofstream outf(fname);
	
	for (int i = 0; i < atomsCount; i++) {
		for (int j = 0; j < signalSize; j++) {
			outf << std::to_string(atoms[i * signalSize + j]) << std::endl;
			if(i==780)
				std::cout << std::to_string(atoms[i * signalSize + j]) << ' ';
		}
		//std::cout << std::endl;
	}
	delete [] atoms;
	return 0;
}
