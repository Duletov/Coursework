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

	auto fname = inputFilePath + ".txt";
	std::ofstream outf(fname);

    for (int i = 0; i < a.getNumSamplesPerChannel(); i++){
        for (int channel = 0; channel < a.getNumChannels(); channel++){
            outf << std::to_string(a.samples[channel][i]) << std::endl;
        }
    }
    std::cout << a.getNumSamplesPerChannel() * a.getNumChannels() << std::endl;
    
    for (int i = 0; i < a.getNumSamplesPerChannel(); i+=10){
        for (int channel = 0; channel < a.getNumChannels(); channel++){
            outf << std::to_string(a.samples[channel][i]) << std::endl;
        }
    }
    std::cout << a.getNumSamplesPerChannel() * a.getNumChannels() / 10;
    
    return 0;
}
