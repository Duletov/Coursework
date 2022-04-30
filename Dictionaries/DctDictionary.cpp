#define _USE_MATH_DEFINES
#include <cmath>

#include "Dictionary.cpp"

class DctDictionary : public Dictionary {
	
	public:
		DctDictionary(int atomsCount, int signalSize, int testSize, int rightBorder) : Dictionary(atomsCount, signalSize, testSize, rightBorder) {}
		
		//Методика построения словаря DCT взята из работы
		//Ahmed, Nasir; Natarajan, T.; Rao, K. R., "Discrete Cosine Transform"
		//IEEE Transactions on Computers, C-23 (1): 90Ц93, January 1974
		void CreateDictionary(double* atoms, double* tests) override {
			
			for (int i = 0; i < atomsCount; i++) {
				for (int j = 0; j < signalSize; j++) {
					atoms[i * signalSize + j] = cos( double(rightBorder) * M_PI / signalSize * (j + 0.5) * i);
				}
			}
			
			for (int i = 0; i < atomsCount; i++) {
				for (int j = 0; j < testSize; j++) {
					tests[i * testSize + j] = cos( double(rightBorder) * M_PI / testSize * (j + 0.5) * i);
				}
			}
		}
};
