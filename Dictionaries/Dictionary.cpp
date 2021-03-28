#pragma once

class Dictionary {
	
	protected:
		int atomsCount, signalSize;
		
		Dictionary(int atomsCount, int signalSize) 
			: atomsCount(atomsCount), signalSize(signalSize) {}	
			
		
	public:		
    	virtual void CreateDictionary(double* atoms) = 0;
};