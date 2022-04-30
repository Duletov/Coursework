#include <iostream>
#include "AudioFile.h"


int main()
{
    // 1. Set a file path to an audio file on your machine
    std::string inputFilePath = "minor.wav";
    
    // 2. Create an AudioFile object and load the audio file
    AudioFile<float> a;
    bool loadedOK = a.load (inputFilePath);
    std::cout << a.getNumSamplesPerChannel() << " ";
    std::cout << a.getNumChannels();
    
    /** If you hit this assert then the file path above
     probably doesn't refer to a valid audio file */
    assert (loadedOK);

	auto fname = inputFilePath + "_1000.txt";
	std::ofstream outf1(fname);
	int k = 0;
	
    for (int i = 0; i < 20000; i++){
        outf1 << std::to_string(a.samples[0][i]) << std::endl;
        k++;
    }
	fname = inputFilePath + "_2000.txt";
	std::ofstream outf2(fname);
    for (int i = 20000; i < 40000; i++){
        outf2 << std::to_string(a.samples[0][i]) << std::endl;
        k++;
    }
    fname = inputFilePath + "_3000.txt";
	std::ofstream outf3(fname);
    for (int i = 40000; i < 60000; i++){
        outf3 << std::to_string(a.samples[0][i]) << std::endl;
        k++;
    }
    fname = inputFilePath + "_4000.txt";
	std::ofstream outf4(fname);
    for (int i = 0; i < 20000; i++){
        outf4 << std::to_string(a.samples[1][i]) << std::endl;
        k++;
    }
    fname = inputFilePath + "_5000.txt";
	std::ofstream outf5(fname);
    for (int i = 20000; i < 40000; i++){
        outf5 << std::to_string(a.samples[1][i]) << std::endl;
        k++;
    }
    fname = inputFilePath + "_6000.txt";
	std::ofstream outf6(fname);
    for (int i = 40000; i < 60000; i++){
        outf6 << std::to_string(a.samples[1][i]) << std::endl;
        k++;
    }
    //std::cout << k << std::endl;
    
    fname = inputFilePath + "_1000_test.txt";
	std::ofstream outf(fname);
	
    for (int i = 0; i < 10000; i++){
        outf << std::to_string(a.samples[0][i]) << std::endl;
    }
    
    return 0;
}
