#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#define sqr(x) ((x)*(x))

#include "Dictionary.cpp"

class MaxTrigSplineDictionary : public Dictionary {
	
	public:
		MaxTrigSplineDictionary(int atomsCount, int signalSize, int testSize, int rightBorder) : Dictionary(atomsCount, signalSize, testSize, rightBorder) {}

		void CreateDictionary(double* atoms, double* tests) override {
			double delta_x = 1.0/(atomsCount - 1.0);
			double xi[atomsCount+4];
			xi[0] = -delta_x;
			for(int i=0;i<atomsCount+3;i++){
				//std::cout << xi[i] << ' ';
				xi[i+1]=xi[i]+delta_x;
			}
			std::cout << xi[atomsCount+3] << std::endl << std::endl;
			for(int i=0;i<atomsCount;i++){
				for(int j=0;j<signalSize;j++){
					double cur_position = 1.0/(signalSize-1)*j+0.5*delta_x;
					//std::cout << cur_position << ' ';
					if((cur_position >= xi[i]) and (cur_position < xi[i+1])){
						atoms[i*signalSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((cur_position - xi[i])/2)))/((sin((xi[i+2]-xi[i])/2))*(sin((xi[i+1]-xi[i])/2)));
						//std::cout << 0 << ' ' << 1 << ' ' << cur_position << ' ' << xi[i] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+1] << std::endl;
					}
					else if((cur_position >= xi[i+1]) and (cur_position < xi[i+2])){
						atoms[i*signalSize+j] = ((cos((xi[i+2]-xi[i+1])/2))/(sin((xi[i+1]-xi[i])/2)))*(((sqr(sin((cur_position-xi[i])/2)))/(sin((xi[i+2]-xi[i])/2)))-((sin((xi[i+3]-xi[i])/2))*(sqr(sin((cur_position-xi[i+1])/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+2]-xi[i+1])/2)))));
						//std::cout << 1 << ' ' << 2 << ' ' << cur_position << ' ' << xi[i+1] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+2] << std::endl;
					}
					else if((cur_position >= xi[i+2]) and (cur_position < xi[i+3])){
						atoms[i*signalSize+j] = (cos((xi[i+2]-xi[i+1])/2))*(sqr(sin((xi[i+3] - cur_position)/2)))/((sin((xi[i+3]-xi[i+1])/2))*(sin((xi[i+3]-xi[i+2])/2)));
						//std::cout << 2 << ' ' << 3 << ' ' << cur_position << ' ' << xi[i+2] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+3] << std::endl;
					}
					else{
						atoms[i*signalSize+j]=0.0;
					}
					//std::cout << atoms[i*signalSize+j] << ' ';
				}
				//std::cout << std::endl;
			}
			return;
		}
};
