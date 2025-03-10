#include <iostream>
#include <vector>
#include <cmath>
#include <complex>

class AudioFingerPrinter {

    private:
        int sampleRate;
        int frameSize;
        int hopSize;
        int numBands;

        std::vector<std::pair<float, float >> bandRanges; 


    public:
        AudioFingerPrinter( int sampleRate = 44100, int frameSize = 2048, int hopSize = 512, int numBands = 32)
        : sampleRate(sampleRate), frameSize(frameSize), hopSize(hopSize), numBands(numBands)
        {

            setupFrequencyBands(300.0, 2000,0); 
        } 

        void setupFrequencyBands(float minFreq, float maxFreq)
        {
            bandRanges.clear(); 
            float logMin = std::log(minFreq); 
            float logMax = std::log(maxFreq); 
            float step = (logMax- logMin) / numBands; 

            for (int i =0; i < numBands; i++)
            {
                float low = std::exp(logMin + i * step); 
                float high = std::exp(logMax + (i +1) * step);
                bandRanges.push_back(std::make_pair(low, high)); 

            }
        }

std::vector<std::vector<float>> computeSpectogram(const std::vector<float> & audioData)
{
    int numFrames = (audioData.size() - frameSize) / hopSize + 1; 
    std::vector<std::vector<float>> spectogram(numFrames, std::vector<float>(frameSize/ 2 + 1, 0.0)); 
    
}
        


}