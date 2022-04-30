#include <fstream>
#include <iostream>

#include "Dictionary.cpp"
#include "GaborDictionary.cpp"
#include "DctDictionary.cpp"
#include "SplineDictionary.cpp"
#include "HypSplineDictionary.cpp"
#include "TrigSplineDictionary.cpp"
#include "MaxTrigSplineDictionary.cpp"
#include "MaxSplineDictionary.cpp"
#include "MaxHypSplineDictionary.cpp"

int main(int argc, char** argv) {
	
	if (argc != 6)	{
		return -1;
	}
	
	int atomsCount =  std::stoi(argv[1]);
	int signalSize = std::stoi(argv[2]);
	int testSize = std::stoi(argv[3]);
	int rightBorder = std::stoi(argv[4]);
	char dictType = argv[5][0];
	
	int temp=0;
	for(int i=4; i<atomsCount;i=i*2)
		temp+=i;
	
	int allAtomsCount = atomsCount + temp;
	
	
	double *atoms = new double[allAtomsCount * signalSize];
	double *tests = new double[allAtomsCount * testSize];
	
	
	switch(dictType){
		case 'd':{
			DctDictionary generator(atomsCount, signalSize, testSize, rightBorder);
			generator.CreateDictionary(atoms, tests);
			break;
		}
		case 'g':{
			GaborDictionary generator(atomsCount, signalSize, testSize, rightBorder);
			generator.CreateDictionary(atoms, tests);
			break;
		}
		case 's':{
			SplineDictionary generator(atomsCount, signalSize, testSize, rightBorder);
			generator.CreateDictionary(atoms, tests);
			break;
		}
		case 'h':{
			HypSplineDictionary generator(atomsCount, signalSize, testSize, rightBorder);
			generator.CreateDictionary(atoms, tests);
			break;
		}
		case 't':{
			TrigSplineDictionary generator(atomsCount, signalSize, testSize, rightBorder);
			generator.CreateDictionary(atoms, tests);
			break;
		}
		case 'a':{
			MaxSplineDictionary generator(atomsCount, signalSize, testSize, rightBorder);
			generator.CreateDictionary(atoms,tests);
			break;
		}
		case 'b':{
			MaxTrigSplineDictionary generator(atomsCount, signalSize, testSize, rightBorder);
			generator.CreateDictionary(atoms,tests);
			break;
		}
		case 'c':{
			MaxHypSplineDictionary generator(atomsCount, signalSize, testSize, rightBorder);
			generator.CreateDictionary(atoms,tests);
			break;
		}
		default:{
			std::cout << "dictType is not defined";
			return 1;
		}
	}
	auto fname = dictType + std::to_string(allAtomsCount) + std::to_string(signalSize) + ".txt";
	std::ofstream outf(fname);
	
	for (int i = 0; i < allAtomsCount; i++) {
		for (int j = 0; j < signalSize; j++) {
			outf << std::to_string(atoms[i * signalSize + j]) << std::endl;
			//std::cout << std::to_string(atoms[i * signalSize + j]) << " ";
		}
	}
	
	if(dictType == 's' or dictType == 't' or dictType == 'h'){
		auto fnamee = dictType + std::to_string(allAtomsCount) + std::to_string(testSize) + ".txt";
		std::ofstream outf1(fnamee);
		
		for (int i = 0; i < allAtomsCount; i++) {
			for (int j = 0; j < testSize; j++) {
				outf1 << std::to_string(tests[i * testSize + j]) << std::endl;
			}
			//std::cout << std::endl;
		}
	}
	
	
	delete [] atoms;
	delete [] tests;
	return 0;
}
