#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <numeric>

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

double normalize(double *atom, int szSignal)
{   
    double normAtom = sqrt(real_inner_product(atom,atom,szSignal));
    
    for ( int i = 0; i < szSignal; i++)
    {
        atom[i] = atom[i]/normAtom;
    }
    return normAtom;
}


class Dictionary{
	protected:
		int nAtoms, szSignal;
		
		Dictionary(int N, int M)
		{
			nAtoms = N;
			szSignal = M;
		}
	public:
    	virtual void CreateDictionary(double *atom) =0;
};


class GaborDictionary: public Dictionary{
	public:
		GaborDictionary(int nAtoms, int szSignal):
			Dictionary(nAtoms, szSignal)
		{
		}
		
		void CreateDictionary(double *atom) override{
			int scale, timeShift;
			double frequency, phase, normAtom, n, m;
			n = double(nAtoms);
			m = double(szSignal);
			scale = 2;
			timeShift = 0;
			frequency = 0.01;
			phase = 0;
			
			for(int i = 0; i<nAtoms; i++){
				for(int j = 0; j<szSignal; j++){
					atom[i*szSignal + j] = (1/sqrt(scale)) * exp(-M_PI*((j/m-i/n)*(j/m-i/n))/(scale*scale)) * cos(2*M_PI*frequency*(j/szSignal-i/nAtoms)+phase);
				}
			}
			//for(int i=0;i<nAtoms;i++){
			//	double normKthAtom = normalize(&atom[i*szSignal],szSignal);
			//}
		}
};

class DCTDictionary: public Dictionary{
	public:
		DCTDictionary(int nAtoms, int szSignal):
			Dictionary(nAtoms, szSignal)
		{
		}
		
		void CreateDictionary(double *atom) override{
			double n = double(nAtoms);
			double m = double(szSignal);
			for (int i=0;i<nAtoms;i++){
				for (int j=0;j<szSignal;j++){
					atom[i*szSignal + j] = cos((M_PI / szSignal) * (j + (.5/szSignal)) * i/nAtoms);
				}
			}
		}
};

int main(int argc, char *argv[]){
	int nAtoms, szSignal;
	nAtoms = stoi(argv[1]);
	szSignal = stoi(argv[2]);
	double atom[nAtoms*szSignal];
	DCTDictionary generator(nAtoms, szSignal);
	generator.CreateDictionary(&atom[0]);
	string fname;
	fname = "d" + to_string(nAtoms) + to_string(szSignal) + ".txt";
	ofstream outf(fname);
	for(int i = 0; i<nAtoms; i++){
		for(int j = 0; j<szSignal; j++){
			outf << to_string(atom[i*szSignal + j]) << endl;
			cout << to_string(atom[i*szSignal + j]) << endl;
		}
	}
}
