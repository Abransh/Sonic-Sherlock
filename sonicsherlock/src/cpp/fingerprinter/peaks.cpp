#include "peaks.h"
#include <algorithm>
#include <cmath>
#include <future>
#include <thread>
#include <unordered_set>

namespace audiofp {

std::vector<Peak> PeakFinder::find_peaks(
    const std::vector<std::vector<double>>& spectrogram,
    int neighborhood_size,
    double threshold,
    int max_peaks) {
    
    std::vector<Peak> peaks;
    
    // Check if spectrogram is empty
    if (spectrogram.empty() || spectrogram[0].empty()) {
        return peaks;
    }
    
    const size_t time_bins = spectrogram.size();
    const size_t freq_bins = spectrogram[0].size();
    
    // Find the maximum value in the spectrogram to determine the threshold
    double max_magnitude = 0.0;
    for (const auto& time_slice : spectrogram) {
        for (double magnitude : time_slice) {
            max_magnitude = std::max(max_magnitude, magnitude);
        }
    }
    
    // Calculate absolute threshold value
    double abs_threshold = max_magnitude * threshold;
    
    // Radius for neighborhood check
    int radius = neighborhood_size / 2;
    if (radius < 1) radius = 1;
    
    // Use parallel processing for large spectrograms
    std::vector<std::future<std::vector<Peak>>> futures;
    const size_t min_chunk_size = 50;  // Minimum time bins per thread
    const size_t max_threads = std::thread::hardware_concurrency();
    const size_t actual_threads = std::min(max_threads, 
                                          (time_bins + min_chunk_size - 1) / min_chunk_size);
    
    if (actual_threads > 1 && time_bins > min_chunk_size) {
        // Parallel peak finding
        const size_t chunk_size = (time_bins + actual_threads - 1) / actual_threads;
        
        for (size_t thread_idx = 0; thread_idx < actual_threads; ++thread_idx) {
            size_t start_time = thread_idx * chunk_size;
            size_t end_time = std::min(start_time + chunk_size, time_bins);
            
            // Skip boundary points that can't be fully processed
            if (thread_idx > 0) {
                start_time += radius;
            }
            if (thread_idx < actual_threads - 1) {
                end_time -= radius;
            }
            
            // Special case for last chunk
            if (thread_idx == actual_threads - 1) {
                end_time = time_bins;
            }
            
            // Lambda to find peaks in a chunk of the spectrogram
            auto find_peaks_in_chunk = [&](size_t start, size_t end) {
                std::vector<Peak> chunk_peaks;
                
                // Skip boundaries for first and last chunks
                size_t effective_start = (start == 0) ? radius : start;
                size_t effective_end = (end == time_bins) ? time_bins - radius : end;
                
                for (size_t t = effective_start; t < effective_end; ++t) {
                    for (size_t f = radius; f < freq_bins - radius; ++f) {
                        double current_value = spectrogram[t][f];
                        
                        // Skip if below threshold
                        if (current_value < abs_threshold) continue;
                        
                        // Check if it's a local maximum in its neighborhood
                        bool is_local_max = true;
                        
                        for (int dt = -radius; dt <= radius && is_local_max; ++dt) {
                            for (int df = -radius; df <= radius; ++df) {
                                // Skip self-comparison
                                if (dt == 0 && df == 0) continue;
                                
                                // Check boundaries
                                size_t nt = t + dt;
                                size_t nf = f + df;
                                
                                if (nt < time_bins && nf < freq_bins && 
                                    spectrogram[nt][nf] >= current_value) {
                                    is_local_max = false;
                                    break;
                                }
                            }
                        }
                        
                        if (is_local_max) {
                            chunk_peaks.emplace_back(t, f, current_value);
                        }
                    }
                }
                
                return chunk_peaks;
            };
            
            futures.push_back(std::async(std::launch::async, 
                                        find_peaks_in_chunk, 
                                        start_time, 
                                        end_time));
        }
        
        // Collect results from all threads
        for (auto& future : futures) {
            auto chunk_peaks = future.get();
            peaks.insert(peaks.end(), chunk_peaks.begin(), chunk_peaks.end());
        }
    } else {
        // Sequential peak finding for smaller spectrograms
        for (size_t t = radius; t < time_bins - radius; ++t) {
            for (size_t f = radius; f < freq_bins - radius; ++f) {
                double current_value = spectrogram[t][f];
                
                // Skip if below threshold
                if (current_value < abs_threshold) continue;
                
                // Check if it's a local maximum in its neighborhood
                bool is_local_max = true;
                
                for (int dt = -radius; dt <= radius && is_local_max; ++dt) {
                    for (int df = -radius; df <= radius; ++df) {
                        // Skip self-comparison
                        if (dt == 0 && df == 0) continue;
                        
                        // Check if any neighbor has a higher value
                        if (spectrogram[t + dt][f + df] >= current_value) {
                            is_local_max = false;
                            break;
                        }
                    }
                }
                
                if (is_local_max) {
                    peaks.emplace_back(t, f, current_value);
                }
            }
        }
    }
    
    // Sort peaks by magnitude (descending)
    std::sort(peaks.begin(), peaks.end());
    
    // Limit the number of peaks if specified
    if (max_peaks > 0 && peaks.size() > static_cast<size_t>(max_peaks)) {
        peaks.resize(max_peaks);
    }
    
    return peaks;
}

std::vector<std::pair<Peak, Peak>> PeakFinder::generate_constellation(
    const std::vector<Peak>& peaks,
    int max_pairs) {
    
    std::vector<std::pair<Peak, Peak>> constellation;
    
    // Constants for constellation map generation
    const int MAX_TIME_DELTA = 100;  // Maximum time difference between peaks
    const int MIN_FREQ_DELTA = 10;   // Minimum frequency difference to avoid duplicate fingerprints
    
    // Create pairs of peaks
    for (size_t i = 0; i < peaks.size(); ++i) {
        const Peak& anchor = peaks[i];
        
        // Try to pair with subsequent peaks within time window
        for (size_t j = i + 1; j < peaks.size(); ++j) {
            const Peak& target = peaks[j];
            
            // Check if target is within time window
            int time_delta = target.time - anchor.time;
            if (time_delta > MAX_TIME_DELTA) {
                // Skip to next anchor as targets are sorted by time
                break;
            }
            
            // Check frequency difference to ensure diversity
            int freq_delta = std::abs(target.frequency - anchor.frequency);
            if (freq_delta < MIN_FREQ_DELTA) {
                continue;
            }
            
            // Add valid pair to constellation map
            constellation.emplace_back(anchor, target);
            
            // Break if we've reached the maximum number of pairs for this anchor
            if (max_pairs > 0 && constellation.size() >= static_cast<size_t>(max_pairs)) {
                break;
            }
        }
    }
    
    return constellation;
}

} // namespace audiofp