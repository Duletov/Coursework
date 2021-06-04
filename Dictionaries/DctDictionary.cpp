#define _USE_MATH_DEFINES
#include <cmath>

#include "Dictionary.cpp"

class DctDictionary : public Dictionary {
	
	public:
		DctDictionary(int atomsCount, int signalSize) : Dictionary(atomsCount, signalSize) {}
		
		//Методика построения словаря DCT взята из работы
		//Ahmed, Nasir; Natarajan, T.; Rao, K. R., "Discrete Cosine Transform"
		//IEEE Transactions on Computers, C-23 (1): 90Ц93, January 1974
		void CreateDictionary(double* atoms) override {
			
			for (int i = 0; i < atomsCount; i++) {
				for (int j = 0; j < signalSize; j++) {
					atoms[i * signalSize + j] = cos( M_PI / signalSize * ((j + 0.5) / signalSize) * i / atomsCount);
				}
			}
		}
};
