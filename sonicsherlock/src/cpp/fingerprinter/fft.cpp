#include "fft.h"
#include <algorithm>
#include <stdexcept>
#include <thread>
#include <future>

namespace audiofp {

std::vector<std::complex<double>> FFT::compute(const std::vector<double>& signal) {
    // Check if signal length is a power of 2
    size_t n = signal.size();
    if (n & (n - 1)) {
        throw std::invalid_argument("Signal length must be a power of 2");
    }
    
    // Convert real signal to complex
    std::vector<std::complex<double>> complex_signal(n);
    for (size_t i = 0; i < n; ++i) {
        complex_signal[i] = std::complex<double>(signal[i], 0.0);
    }
    
    // Use iterative FFT for better performance
    const double PI = 3.14159265358979323846;
    
    // Bit reversal
    bit_reverse(complex_signal);
    
    // Cooley-Tukey FFT algorithm (iterative implementation)
    for (size_t s = 1; s < n; s *= 2) {
        size_t m = s * 2;
        std::complex<double> wm = std::polar(1.0, -2.0 * PI / m);
        
        // Process each butterfly group in parallel if large enough
        if (n >= 4096 && s >= 64) {
            size_t num_threads = std::min(std::thread::hardware_concurrency(), n / m);
            std::vector<std::future<void>> futures;
            
            for (size_t t = 0; t < num_threads; ++t) {
                size_t start = t * (n / num_threads);
                size_t end = (t == num_threads - 1) ? n : (t + 1) * (n / num_threads);
                start = (start / m) * m; // Align to butterfly boundaries
                
                futures.push_back(std::async(std::launch::async, [&, start, end]() {
                    for (size_t k = start; k < end; k += m) {
                        std::complex<double> w(1.0, 0.0);
                        for (size_t j = 0; j < s; ++j) {
                            std::complex<double> t = w * complex_signal[k + j + s];
                            std::complex<double> u = complex_signal[k + j];
                            complex_signal[k + j] = u + t;
                            complex_signal[k + j + s] = u - t;
                            w *= wm;
                        }
                    }
                }));
            }
            
            for (auto& future : futures) {
                future.wait();
            }
        } else {
            // Sequential processing for smaller transforms
            for (size_t k = 0; k < n; k += m) {
                std::complex<double> w(1.0, 0.0);
                for (size_t j = 0; j < s; ++j) {
                    std::complex<double> t = w * complex_signal[k + j + s];
                    std::complex<double> u = complex_signal[k + j];
                    complex_signal[k + j] = u + t;
                    complex_signal[k + j + s] = u - t;
                    w *= wm;
                }
            }
        }
    }
    
    return complex_signal;
}

std::vector<double> FFT::magnitude_spectrum(const std::vector<std::complex<double>>& fft_result) {
    std::vector<double> magnitudes(fft_result.size());
    
    // Use parallel processing for large arrays
    if (fft_result.size() >= 4096) {
        size_t num_threads = std::min(std::thread::hardware_concurrency(), fft_result.size() / 1024);
        std::vector<std::future<void>> futures;
        
        for (size_t t = 0; t < num_threads; ++t) {
            size_t start = t * (fft_result.size() / num_threads);
            size_t end = (t == num_threads - 1) ? fft_result.size() : (t + 1) * (fft_result.size() / num_threads);
            
            futures.push_back(std::async(std::launch::async, [&, start, end]() {
                for (size_t i = start; i < end; ++i) {
                    // Optimized magnitude calculation
                    double real = fft_result[i].real();
                    double imag = fft_result[i].imag();
                    magnitudes[i] = std::sqrt(real * real + imag * imag);
                }
            }));
        }
        
        for (auto& future : futures) {
            future.wait();
        }
    } else {
        // Sequential processing for smaller arrays
        for (size_t i = 0; i < fft_result.size(); ++i) {
            double real = fft_result[i].real();
            double imag = fft_result[i].imag();
            magnitudes[i] = std::sqrt(real * real + imag * imag);
        }
    }
    
    return magnitudes;
}

std::vector<double> FFT::apply_window(const std::vector<double>& signal, int window_type) {
    std::vector<double> windowed_signal(signal.size());
    const double PI = 3.14159265358979323846;
    
    switch (window_type) {
        case 1: // Hanning window
            for (size_t i = 0; i < signal.size(); ++i) {
                double multiplier = 0.5 * (1.0 - std::cos(2.0 * PI * i / (signal.size() - 1)));
                windowed_signal[i] = signal[i] * multiplier;
            }
            break;
        case 2: // Hamming window
            for (size_t i = 0; i < signal.size(); ++i) {
                double multiplier = 0.54 - 0.46 * std::cos(2.0 * PI * i / (signal.size() - 1));
                windowed_signal[i] = signal[i] * multiplier;
            }
            break;
        default: // Rectangular window (no change)
            windowed_signal = signal;
            break;
    }
    
    return windowed_signal;
}

void FFT::bit_reverse(std::vector<std::complex<double>>& signal) {
    size_t n = signal.size();
    size_t j = 0;
    
    for (size_t i = 0; i < n - 1; ++i) {
        if (i < j) {
            std::swap(signal[i], signal[j]);
        }
        
        size_t m = n >> 1;
        while (j >= m && m > 0) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
}

} // namespace audiofp