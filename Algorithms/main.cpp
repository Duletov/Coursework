#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>
#include <cmath>
#include <boost/spirit/include/qi_real.hpp>

#include "Algorithm.cpp"
#include "OMP.cpp"
#include "OMP_str.cpp"
#include "Lasso.cpp"
#include "StOMP.cpp"
#include "SOMP.cpp"

int main(int argc, char** argv) {
	//std::ios_base::sync_with_stdio(false);
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
	
	
	fname = "signals/Black_Hole_Billiards.wav_boost.txt";
	std::ifstream in2(fname.c_str());
	for(int i=0;i<szSignal;i++){
		vSignal[i] = atan(10*i*1.0/szSignal);
		//in2 >> vSignal[i];
	}
	fname = "signals/Black_Hole_Billiards.wav_boost_test.txt";
	std::ifstream in3(fname.c_str());
	for(int i=0;i<szTest;i++){
		rSignal[i] = atan(10*i*1.0/szTest);
		//in3 >> rSignal[i];
	}
	std::cout << std::endl;
	auto start = clock();
	
	switch(algoType){
		case 'o':{
			OMPAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		case 'l':{
			LassoAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		case 's':{
			StOMPAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		case 't':{
			OMPsAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		case 'p':{
			SplineOMPAlgorithm instance(nAtoms, szSignal, szTest);
			instance.RunAlgorithm(vSignal, rSignal, mDictionary, fullDictionary);
			break;
		}
		default:{
			std::cout << "dictType is not defined";
			return 1;
		}
	}
	auto end = clock();
	
	std::cout << (end - start) / CLOCKS_PER_SEC;
	
		    std::cout << "ok";
	return 0;
}
