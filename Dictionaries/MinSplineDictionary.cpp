#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#define sqr(x) ((x)*(x))

#include "Dictionary.cpp"

class MinSplineDictionary : public Dictionary {
	
	public:
		MinSplineDictionary(int atomsCount, int signalSize) : Dictionary(atomsCount, signalSize) {}

		void CreateDictionary(double* atoms) override {
			double delta_x = 1.0/(atomsCount - 1.0);
			double xi[atomsCount+4];
			xi[0] = -1.5*delta_x;
			for(int i=0;i<atomsCount+3;i++){
				//std::cout << xi[i] << ' ';
				xi[i+1]=xi[i]+delta_x;
			}
			std::cout << xi[atomsCount+3] << std::endl << std::endl;
			for(int i=0;i<atomsCount;i++){
				for(int j=0;j<signalSize;j++){
					double cur_position = 1.0/(signalSize-1)*j;
					//std::cout << cur_position << ' ';
					if((cur_position >= xi[i]) and (cur_position < xi[i+1])){
						atoms[i*signalSize+j] = ((cosh((xi[i+2]-xi[i+1])/2))*sqr(sinh((cur_position-xi[i])/2)))/(sinh((xi[i+1]-xi[i])/2)*sinh((xi[i+2]-xi[i])/2));
						//std::cout << 0 << ' ' << 1 << ' ' << cur_position << ' ' << xi[i] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+1] << std::endl;
					}
					else if((cur_position >= xi[i+1]) and (cur_position < xi[i+2])){
						atoms[i*signalSize+j] = ((cosh((xi[i+2]-xi[i+1])/2))/(sinh((xi[i+1]-xi[i])/2)))*(((sqr(sinh((cur_position-xi[i])/2)))/(sinh((xi[i+2]-xi[i])/2)))-(((sinh((xi[i+3]-xi[i])/2))*(sqr(sinh((cur_position-xi[i+1])/2))))/((sinh((xi[i+3]-xi[i+1])/2))*(sinh((xi[i+2]-xi[i+1])/2)))));
						//std::cout << 1 << ' ' << 2 << ' ' << cur_position << ' ' << xi[i+1] << ' ' << atoms[i*signalSize+j] << ' ' << xi[i+2] << std::endl;
					}
					else if((cur_position >= xi[i+2]) and (cur_position < xi[i+3])){
						atoms[i*signalSize+j] = ((cosh((xi[i+2]-xi[i+1])/2))*(sqr(sinh((xi[i+3]-cur_position)/2))))/((sinh((xi[i+3]-xi[i+1])/2))*(sinh((xi[i+3]-xi[i+2])/2)));
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
