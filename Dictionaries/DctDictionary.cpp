#pragma once
#include <cmath>

class DctDictionary : public Dictionary {
	
	public:
		DctDictionary(int atomsCount, int signalSize) : Dictionary(atomsCount, signalSize) {}
		
		//Методика построения DCT словаря взята ...
		void CreateDictionary(double* atoms) override {
			
			for (int i = 0; i < atomsCount; i++) {
				for (int j = 0; j < signalSize; j++) {
					atoms[i * signalSize + j] = cos( M_PI / signalSize * (j + (0.5 / signalSize)) * i / atomsCount);
				}
			}
		}
};