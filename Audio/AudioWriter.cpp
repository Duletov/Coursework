#include "AudioFile.h"
#include <fstream>


int main(){
	AudioFile<double> a;
	
	a.setNumChannels (2);
    a.setNumSamplesPerChannel (60000);
	
	const float sampleRate = 44100.f;
    const float frequencyInHz = 440.f;
	
	std::string fname = "imggg1.txt";
	std::ifstream in0(fname.c_str());
	for(int i=0;i<60000;i++){
		in0 >> a.samples[0][i];
	}
	
	fname = "imggg1.txt";
	std::ifstream in1(fname.c_str());
	for(int i=0;i<60000;i++){
		in1 >> a.samples[1][i];
	}
	
	std::string filePath = "sine-wave.wav"; // change this to somewhere useful for you
    a.save ("sine-wave.wav", AudioFileFormat::Wave);
}
