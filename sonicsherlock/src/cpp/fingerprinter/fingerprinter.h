#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include "peaks.h"

namespace audiofp {

// Hash structure for fingerprints
struct FingerprintHash {
    uint32_t hash;         // 32-bit hash value
    uint32_t time_offset;  // Time offset for matching
    
    FingerprintHash(uint32_t h, uint32_t t) : hash(h), time_offset(t) {}
};

class Fingerprinter {
public:
    /**
     * Generate a spectrogram from audio samples
     * 
     * @param samples Audio samples (mono)
     * @param sample_rate Sample rate in Hz
     * @param window_size FFT window size
     * @param hop_size Hop size between consecutive FFT windows
     * @return 2D spectrogram (time x frequency)
     */
    static std::vector<std::vector<double>> generate_spectrogram(
        const std::vector<double>& samples,
        int sample_rate,
        int window_size = 1024,
        int hop_size = 512
    );
    
    /**
     * Generate audio fingerprints from a spectrogram
     * 
     * @param spectrogram 2D spectrogram
     * @param sample_rate Original audio sample rate (for scaling)
     * @return Vector of hash fingerprints
     */
    static std::vector<FingerprintHash> generate_fingerprints(
        const std::vector<std::vector<double>>& spectrogram,
        int sample_rate
    );
    
    /**
     * Generate fingerprints directly from audio samples
     * 
     * @param samples Audio samples (mono)
     * @param sample_rate Sample rate in Hz
     * @return Vector of hash fingerprints
     */
    static std::vector<FingerprintHash> fingerprint_audio(
        const std::vector<double>& samples,
        int sample_rate
    );
    
    /**
     * Match fingerprints against a database (to be implemented in Python)
     * This is just a C++ interface for Python binding
     * 
     * @param fingerprints Query fingerprints
     * @return Match confidence (0.0-1.0)
     */
    static double match_fingerprints(
        const std::vector<FingerprintHash>& fingerprints,
        const std::vector<FingerprintHash>& db_fingerprints
    );

private:
    // Constants for fingerprinting
    static constexpr int TARGET_ZONE_SIZE = 10;  // Number of points in target zone
    static constexpr int FREQUENCY_BANDS = 6;    // Number of frequency bands
    
    // Hash function for peak pairs
    static uint32_t create_hash(
        const Peak& anchor, 
        const Peak& target,
        int frequency_bands
    );
};

} // namespace audiofp