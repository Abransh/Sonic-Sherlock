#pragma once

#include <vector>
#include <complex>
#include <cmath>

namespace audiofp {

class FFT {
public:
    /**
     * Compute the Fast Fourier Transform of the input signal
     * 
     * @param signal Input time-domain signal
     * @return Complex-valued frequency-domain representation
     */
    static std::vector<std::complex<double>> compute(const std::vector<double>& signal);
    
    /**
     * Computes the magnitude spectrum from the complex FFT result
     * 
     * @param fft_result Complex FFT result
     * @return Magnitude spectrum (absolute values)
     */
    static std::vector<double> magnitude_spectrum(const std::vector<std::complex<double>>& fft_result);
    
    /**
     * Applies a window function to the signal before FFT to reduce spectral leakage
     * 
     * @param signal Input signal to be windowed
     * @param window_type Type of window (0=rectangular, 1=hanning, 2=hamming)
     * @return Windowed signal
     */
    static std::vector<double> apply_window(const std::vector<double>& signal, int window_type);

private:
    // Fast recursive implementation of the Cooley-Tukey FFT algorithm
    static std::vector<std::complex<double>> fft_recursive(const std::vector<std::complex<double>>& signal);
    
    // Bit reversal for FFT optimization
    static void bit_reverse(std::vector<std::complex<double>>& signal);
};

} // namespace audiofp