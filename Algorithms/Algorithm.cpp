#pragma once

class Algorithm{
	
	protected:
		int nAtoms, szSignal, szTest;
		
		Algorithm(int nAtoms, int szSignal, int szTest) 
			: nAtoms(nAtoms), szSignal(szSignal),  szTest(szTest){}	
			
		
	public:		
    	virtual void RunAlgorithm(double* vSignal, double* rSignal, double* mDictionary, double* fullDictionary) = 0;
};
