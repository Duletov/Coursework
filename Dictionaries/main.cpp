#include <fstream>
#include <iostream>

#include "Dictionary.cpp"
#include "GaborDictionary.cpp"
#include "DctDictionary.cpp"
#include "SplineDictionary.cpp"

int main(int argc, char** argv) {
	
	if (argc != 4)	{
		return -1;
	}
	
	int atomsCount =  std::stoi(argv[1]);
	int signalSize = std::stoi(argv[2]);
	char dictType = argv[3][0];
	
	double atoms[atomsCount * signalSize];
	
	if (dictType == 'd'){
		DctDictionary generator(atomsCount, signalSize);
		generator.CreateDictionary(atoms);
	}
	else{
		if (dictType == 'g'){
			GaborDictionary generator(atomsCount, signalSize);
			generator.CreateDictionary(atoms);
		}
		else{
			if (dictType == 's'){
				SplineDictionary generator(atomsCount, signalSize);
				generator.CreateDictionary(atoms);
			}
			else{
				std::cout << "dictType is not defined";
				return 1;
			}
		}
	}
	
	auto fname = dictType + std::to_string(atomsCount) + std::to_string(signalSize) + ".txt";
	std::ofstream outf(fname);
	
	for (int i = 0; i < atomsCount; i++) {
		for (int j = 0; j < signalSize; j++) {
			outf << std::to_string(atoms[i * signalSize + j]) << std::endl;
		}
	}
	
	return 0;
}
