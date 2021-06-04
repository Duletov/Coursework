#define _USE_MATH_DEFINES
#include <cmath>
#define sqr(x) ((x)*(x))

#include "Dictionary.cpp"

class GaborDictionary : public Dictionary {
	
	public:
		GaborDictionary(int atomsCount, int signalSize) : Dictionary(atomsCount, signalSize) {}
		
		//Методика построения словаря Габора, равно как и значения параметров, взяты из работы
		//S. Chu, S. Narayanan, C.-C. Jay Kuo "Environmental Sound Recognition With Time–Frequency Audio Features"  
		//IEEE TRANSACTIONS ON AUDIO, SPEECH, AND LANGUAGE PROCESSING, VOL. 17, NO. 6, AUGUST 2009
		void CreateDictionary(double* atoms) override {
			
			int scale = 2;
			double frequency = 0.01; //откуда это значение?
			double phase = 0;
			
			for(int i = 0; i < atomsCount; i++) {
				for(int j = 0; j < signalSize; j++) {
					
					auto parameter = 1.0 * j / signalSize - 1.0 * i / atomsCount; //откуда это значение?
					auto exponent = exp(- M_PI * sqr(parameter) / sqr(scale) );
					auto cosinus = cos(2 * M_PI * frequency * parameter + phase);
					
					atoms[i * signalSize + j] = 1.0 / sqrt(scale)  * exponent * cosinus;
				}
			}
			
			for(int i = 0; i < atomsCount; i++) {
				//TODO: подключить библиотеку MathUtils и вызвать функцию оттуда
				//normalize_vector(&atoms[i * signalSize], signalSize);
			}
		}
};
