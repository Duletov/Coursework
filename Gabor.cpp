#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

double real_inner_product(double *v1, double *v2, int szVector)
{
  double sum = 0;
  for ( int i = 0; i < szVector; i++)
  {
    sum += v1[i]*v2[i];    
  }
  
  return sum;
}

int main(int argc, char *argv[]){
	int i, j, nAtoms, szSignal, scale, timeShift;
	double frequency, phase, normAtom, tmp, n, m;
	string fname;
	nAtoms = stoi(argv[1]);
	szSignal = stoi(argv[2]);
	n = double(nAtoms);
	m = double(szSignal);
	fname = "g" + to_string(nAtoms) + to_string(szSignal) + ".txt";
	scale = 2;
	timeShift = 0;
	frequency = 0.01;
	phase = 0;
	
	double atom[nAtoms][szSignal];
	for(i = 0; i<nAtoms; i++){
		for(j = 0; j<szSignal; j++){
			atom[i][j] = (1/sqrt(scale)) * exp(-M_PI*((j/m-i/n)*(j/m-i/n))/(scale*scale)) * cos(2*M_PI*frequency*(j/szSignal-i/nAtoms)+phase);
			printf("%f\n", atom[i][j]);
		}
	}
	ofstream outf(fname);
	for(i = 0; i<nAtoms; i++){
		for(j = 0; j<szSignal; j++){
			outf << to_string(atom[i][j]) << endl;
		}
	}
}
