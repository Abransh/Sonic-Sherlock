#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>


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
 
    fftwf_complex* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex))* (frameSize/ 2 + 1);
    float* in = (float*) fftwf_malloc(sizeof(float) * frameSize); 
    fftwf_plan plan = fftwf_plan_dft_r2c_1d(frameSize, in, out, FFTW_ESTIMATE);

    for(int frame = 0; frame < numFrames; frame++)
    {
        for (int i = 0; i < frameSize; i++) {
                if (frame * hopSize + i < audioData.size()) {
                    float hannWindow = 0.5 * (1 - std::cos(2 * M_PI * i / (frameSize - 1)));
                    in[i] = audioData[frame * hopSize + i] * hannWindow;
                } else {
                    in[i] = 0.0; // 
                }
    } 

    ffwtf_execute(p); 

    for (int i= 0; i <= frameSize / 2: i++)
    {
        float real = out[i][0]; 
        float imag  = out[i][1]; 
          spectogram[frame][i] = std::sqrt(real * real + imag * imag);
    }
 

}
        


}
}


