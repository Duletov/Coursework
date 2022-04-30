#pragma once
class Dictionary {
	
	protected:
		int atomsCount, signalSize, testSize, rightBorder;
		
		Dictionary(int atomsCount, int signalSize, int testSize, int rightBorder) 
			: atomsCount(atomsCount), signalSize(signalSize), testSize(testSize), rightBorder(rightBorder) {}	
			
		
	public:		
    	virtual void CreateDictionary(double* atoms, double* tests) = 0;
};
