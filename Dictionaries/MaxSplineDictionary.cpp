#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#define sqr(x) ((x)*(x))

#include "Dictionary.cpp"

class MaxSplineDictionary : public Dictionary {
	
	public:
		MaxSplineDictionary(int atomsCount, int signalSize, int testSize, int rightBorder) : Dictionary(atomsCount, signalSize, testSize, rightBorder) {}

		void CreateDictionary(double* atoms, double* tests) override {
			double delta_x = double(rightBorder)/(atomsCount - 1.0);
			double xi[atomsCount+4];
			xi[0] = -0.5*delta_x;
			xi[1] = 0.0;
			xi[atomsCount+3] = 1.0;
			xi[atomsCount+2] = 1.0 + 0.5*delta_x;
			for(int i=1;i<atomsCount+1;i++){
				//std::cout << xi[i] << ' ';
				xi[i+1]=xi[i]+delta_x;
			}
			//std::cout << xi[atomsCount+3] << std::endl << std::endl;
			for(int i=0;i<atomsCount;i++){
				for(int j=0;j<signalSize;j++){
					double cur_position = (1.0/(signalSize-1)*j)+0.5*delta_x;
					//std::cout << cur_position << ' ';
					if((cur_position >= xi[i]) and (cur_position < xi[i+1])){
						atoms[i*signalSize+j] = sqr(cur_position - xi[i])/((xi[i+1]-xi[i])*(xi[i+2]-xi[i]));
						//std::cout << 0 << ' ' << 1 << ' ' << cur_position << ' ' << xi[i] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+1] << std::endl;
					}
					else if((cur_position >= xi[i+1]) and (cur_position < xi[i+2])){
						atoms[i*signalSize+j] = (1.0/(xi[i+1]-xi[i]))*(((sqr(cur_position - xi[i]))/(xi[i+2]-xi[i]))-(sqr(cur_position - xi[i+1])*(xi[i+3]-xi[i]))/((xi[i+2]-xi[i+1])*(xi[i+3]-xi[i+1])));
						//std::cout << 1 << ' ' << 2 << ' ' << cur_position << ' ' << xi[i+1] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+2] << std::endl;
					}
					else if((cur_position >= xi[i+2]) and (cur_position < xi[i+3])){
						atoms[i*signalSize+j] = sqr(cur_position - xi[i+3])/((xi[i+3]-xi[i+1])*(xi[i+3]-xi[i+2]));
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
