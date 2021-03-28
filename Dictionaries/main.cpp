#include <fstream>

#include "Dictionary.cpp"
#include "GaborDictionary.cpp"
#include "DctDictionary.cpp"

int main(int argc, char** argv) {
	
	if (argc != 4)	{
		return -1;
	}
	
	int atomsCount =  std::stoi(argv[1]);
	int signalSize = std::stoi(argv[2]);
	int dictType = std::stoi(argv[3]);
	
	//TODO: здесь нужно сделать создание объекта нужного типа в зависимости от dictType
	DctDictionary generator(atomsCount, signalSize);
	
	double atoms[atomsCount * signalSize];
	generator.CreateDictionary(&atoms[0]);
	
	auto fname = "d" + std::to_string(atomsCount) + std::to_string(signalSize) + ".txt";
	std::ofstream outf(fname);
	
	for (int i = 0; i < atomsCount; i++) {
		for (int j = 0; j < signalSize; j++) {
			outf << std::to_string(atoms[i * signalSize + j]) << std::endl;
		}
	}
	
	return 0;
}