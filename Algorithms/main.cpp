#include <fstream>
#include <iostream>
#include <string>

#include "Algorithm.cpp"
#include "OMP.cpp"
#include "Lasso.cpp"
#include "StOMP.cpp"

int main(int argc, char** argv) {
	if (argc != 6)	{
		return -1;
	}
	
	int nAtoms = std::stoi(argv[1]);
	int szSignal = std::stoi(argv[2]);
	int szTest = std::stoi(argv[3]);
	char dictType = argv[4][0];
	char algoType = argv[5][0];
	
	double *vSignal = new double[szSignal];
	double *rSignal = new double[szTest];
	double *mDictionary = new double[szSignal*nAtoms];
	double *fullDictionary = new double[nAtoms*szTest];
	
	if(!vSignal || !rSignal || !mDictionary || !fullDictionary)
		exit(2);
	
	std::string fname = "dicts/" + std::string(1, dictType) + std::to_string(nAtoms) + std::to_string(szSignal) + ".txt";
	std::ifstream in(fname.c_str());
	if (in.is_open()){
		for(int i=0;i<szSignal*nAtoms;i++){
			in >> mDictionary[i];
		}
	}
	in.close();
	
	fname = "dicts/" + std::string(1, dictType) + std::to_string(nAtoms) + std::to_string(szTest) + ".txt";
	std::ifstream inn(fname.c_str());
	if (inn.is_open()){
		for(int i=0;i<szTest*nAtoms;i++){
			inn >> fullDictionary[i];
		}
	}
	inn.close();
	
	
	fname = "signals/Black_Hole_Billiards.wav_1000.txt";
	std::ifstream in2(fname.c_str());
	for(int i=0;i<szSignal;i++){
		in2 >> vSignal[i];
	}
	
	fname = "signals/Black_Hole_Billiards.wav_1000_test.txt";
	std::ifstream in3(fname.c_str());
	for(int i=0;i<szTest;i++){
		in3 >> rSignal[i];
	}
	
	
	if (algoType == 'o'){
		OMPAlgorithm instance(nAtoms, szSignal, szTest);
		instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
	}
	else{
		if (algoType == 'l'){
			LassoAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
		}
		else{
			if (algoType == 's'){
				StOMPAlgorithm instance(nAtoms, szSignal, szTest);
				instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			}
			else{
				std::cout << "algoType is not defined";
				return 1;
			}
		}
	}
	
	return 0;
}
