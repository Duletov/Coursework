#include <iostream>
#include "AudioFile.h"


int main()
{
    // 1. Set a file path to an audio file on your machine
    std::string inputFilePath = "Black_Hole_Billiards.wav";
    
    // 2. Create an AudioFile object and load the audio file
    AudioFile<float> a;
    bool loadedOK = a.load (inputFilePath);
    
    /** If you hit this assert then the file path above
     probably doesn't refer to a valid audio file */
    assert (loadedOK);

	auto fname = inputFilePath + "_1000.txt";
	std::ofstream outf1(fname);
	int k = 0;
	
    for (int i = 1; i < 10001; i+=10){
        outf1 << std::to_string(a.samples[0][i]) << std::endl;
        k++;
    }
    std::cout << k << std::endl;
    
    fname = inputFilePath + "_1000_test.txt";
	std::ofstream outf2(fname);
	
    for (int i = 0; i < 10000; i++){
        outf2 << std::to_string(a.samples[0][i]) << std::endl;
    }
    std::cout << a.getNumSamplesPerChannel() * a.getNumChannels() / 100;
    
    return 0;
}
