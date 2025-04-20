#include "fingerprinter.h"
#include "fft.h"
#include "peaks.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <future>
#include <thread>
#include <unordered_map>

namespace audiofp {

std::vector<std::vector<double>> Fingerprinter::generate_spectrogram(
    const std::vector<double>& samples,
    int sample_rate,
    int window_size,
    int hop_size) {
    
    // Check for valid input parameters
    if (samples.empty() || window_size <= 0 || hop_size <= 0) {
        return {};
    }
    
    // Calculate the number of time frames
    const size_t num_samples = samples.size();
    const size_t num_frames = 1 + (num_samples - window_size) / hop_size;
    
    // Check if we have enough samples for at least one frame
    if (num_frames == 0) {
        return {};
    }
    
    // Prepare the spectrogram matrix (time bins x frequency bins)
    // Frequency bins will be window_size/2 + 1 (as we use real FFT)
    const size_t num_freq_bins = window_size / 2 + 1;
    std::vector<std::vector<double>> spectrogram(num_frames, std::vector<double>(num_freq_bins, 0.0));
    
    // Use parallel processing for large spectrograms
    const size_t min_chunk_size = 50;  // Minimum frames per thread
    const size_t max_threads = std::thread::hardware_concurrency();
    const size_t actual_threads = std::min(max_threads, 
                                        (num_frames + min_chunk_size - 1) / min_chunk_size);
    
    if (actual_threads > 1 && num_frames > min_chunk_size) {
        // Parallel spectrogram generation
        const size_t chunk_size = (num_frames + actual_threads - 1) / actual_threads;
        std::vector<std::future<void>> futures;
        
        for (size_t thread_idx = 0; thread_idx < actual_threads; ++thread_idx) {
            size_t start_frame = thread_idx * chunk_size;
            size_t end_frame = std::min(start_frame + chunk_size, num_frames);
            
            // Process a chunk of frames in parallel
            futures.push_back(std::async(std::launch::async, [&](size_t start, size_t end) {
                // Create a temporary buffer for the current window
                std::vector<double> windowed_buffer(window_size);
                
                for (size_t frame = start; frame < end; ++frame) {
                    // Extract the current window from samples
                    size_t offset = frame * hop_size;
                    
                    // Copy samples for the current window
                    for (size_t i = 0; i < window_size && (offset + i) < num_samples; ++i) {
                        windowed_buffer[i] = samples[offset + i];
                    }
                    
                    // Apply a Hamming window to reduce spectral leakage
                    windowed_buffer = FFT::apply_window(windowed_buffer, 2);  // 2 = Hamming window
                    
                    // Compute FFT
                    auto fft_result = FFT::compute(windowed_buffer);
                    
                    // Compute magnitude spectrum
                    auto magnitudes = FFT::magnitude_spectrum(fft_result);
                    
                    // Store the first half of the spectrum (plus DC) in the spectrogram
                    for (size_t i = 0; i < num_freq_bins; ++i) {
                        spectrogram[frame][i] = magnitudes[i];
                    }
                }
            }, start_frame, end_frame));
        }
        
        // Wait for all threads to complete
        for (auto& future : futures) {
            future.wait();
        }
    } else {
        // Sequential processing for smaller spectrograms
        // Create a buffer for the current window
        std::vector<double> windowed_buffer(window_size);
        
        for (size_t frame = 0; frame < num_frames; ++frame) {
            // Extract the current window from samples
            size_t offset = frame * hop_size;
            
            // Copy samples for the current window
            for (size_t i = 0; i < window_size && (offset + i) < num_samples; ++i) {
                windowed_buffer[i] = samples[offset + i];
            }
            
            // Apply a Hamming window to reduce spectral leakage
            windowed_buffer = FFT::apply_window(windowed_buffer, 2);  // 2 = Hamming window
            
            // Compute FFT
            auto fft_result = FFT::compute(windowed_buffer);
            
            // Compute magnitude spectrum
            auto magnitudes = FFT::magnitude_spectrum(fft_result);
            
            // Store the first half of the spectrum (plus DC) in the spectrogram
            for (size_t i = 0; i < num_freq_bins; ++i) {
                spectrogram[frame][i] = magnitudes[i];
            }
        }
    }
    
    return spectrogram;
}

std::vector<FingerprintHash> Fingerprinter::generate_fingerprints(
    const std::vector<std::vector<double>>& spectrogram,
    int sample_rate) {
    
    std::vector<FingerprintHash> fingerprints;
    
    // Check for empty spectrogram
    if (spectrogram.empty() || spectrogram[0].empty()) {
        return fingerprints;
    }
    
    // Find spectral peaks in the spectrogram
    auto peaks = PeakFinder::find_peaks(spectrogram, 10, 0.3, 1000);
    
    // Generate constellation map from peaks
    auto constellation = PeakFinder::generate_constellation(peaks, 0);
    
    // Reserve space for fingerprints
    fingerprints.reserve(constellation.size());
    
    // Create fingerprint hashes from peak pairs
    for (const auto& peak_pair : constellation) {
        const Peak& anchor = peak_pair.first;
        const Peak& target = peak_pair.second;
        
        // Create hash from the peak pair
        uint32_t hash = create_hash(anchor, target, FREQUENCY_BANDS);
        
        // Time offset (in frames) for matching
        uint32_t time_offset = anchor.time;
        
        // Add to fingerprints
        fingerprints.emplace_back(hash, time_offset);
    }
    
    return fingerprints;
}

std::vector<FingerprintHash> Fingerprinter::fingerprint_audio(
    const std::vector<double>& samples,
    int sample_rate) {
    
    // Generate spectrogram from audio samples
    auto spectrogram = generate_spectrogram(samples, sample_rate);
    
    // Generate fingerprints from spectrogram
    return generate_fingerprints(spectrogram, sample_rate);
}

double Fingerprinter::match_fingerprints(
    const std::vector<FingerprintHash>& query_fingerprints,
    const std::vector<FingerprintHash>& db_fingerprints) {
    
    // Simple matching algorithm that counts exact hash matches
    std::unordered_map<uint32_t, uint32_t> query_hash_map;
    
    // Index query fingerprints by hash
    for (const auto& fp : query_fingerprints) {
        query_hash_map[fp.hash] = fp.time_offset;
    }
    
    // Count matches
    size_t match_count = 0;
    std::unordered_map<int, int> time_deltas;
    
    for (const auto& db_fp : db_fingerprints) {
        auto it = query_hash_map.find(db_fp.hash);
        if (it != query_hash_map.end()) {
            // Found a matching hash
            match_count++;
            
            // Calculate time delta for alignment
            int time_delta = static_cast<int>(db_fp.time_offset) - static_cast<int>(it->second);
            time_deltas[time_delta]++;
        }
    }
    
    // Find the most common time delta (best alignment)
    int max_delta_count = 0;
    for (const auto& delta : time_deltas) {
        if (delta.second > max_delta_count) {
            max_delta_count = delta.second;
        }
    }
    
    // Calculate confidence score (0.0 - 1.0)
    double confidence = 0.0;
    if (!query_fingerprints.empty()) {
        confidence = static_cast<double>(max_delta_count) / query_fingerprints.size();
    }
    
    return confidence;
}

uint32_t Fingerprinter::create_hash(
    const Peak& anchor,
    const Peak& target,
    int frequency_bands) {
    
    // Quantize frequency values to reduce hash variations
    uint32_t freq1 = static_cast<uint32_t>(anchor.frequency) % frequency_bands;
    uint32_t freq2 = static_cast<uint32_t>(target.frequency) % frequency_bands;
    
    // Calculate time delta between points
    uint32_t time_delta = static_cast<uint32_t>(target.time - anchor.time);
    
    // Combine values into a 32-bit hash
    // Format: [freq1(3 bits)][freq2(3 bits)][time_delta(26 bits)]
    uint32_t hash = (freq1 << 29) | (freq2 << 26) | (time_delta & 0x3FFFFFF);
    
    return hash;
}

} // namespace audiofp