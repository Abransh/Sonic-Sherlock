#pragma once

#include <vector>
#include <utility>

namespace audiofp {

struct Peak {
    int time;       // Time bin
    int frequency;  // Frequency bin
    double magnitude; // Peak magnitude
    
    Peak(int t, int f, double m) : time(t), frequency(f), magnitude(m) {}
    
    // Sorting operator for peaks (by magnitude, descending)
    bool operator<(const Peak& other) const {
        return magnitude > other.magnitude; // Note: reversed for descending sort
    }
};

class PeakFinder {
public:
    /**
     * Find spectral peaks in a spectrogram
     * 
     * @param spectrogram 2D spectrogram (time x frequency)
     * @param neighborhood_size Size of the local neighborhood to check (default: 3x3)
     * @param threshold Minimum magnitude threshold (relative to max value)
     * @param max_peaks Maximum number of peaks to return (0 = unlimited)
     * @return Vector of peak information (time, frequency, magnitude)
     */
    static std::vector<Peak> find_peaks(
        const std::vector<std::vector<double>>& spectrogram,
        int neighborhood_size = 3,
        double threshold = 0.3,
        int max_peaks = 0
    );
    
    /**
     * Generates a constellation map from peak data
     * 
     * @param peaks Detected spectral peaks
     * @param max_pairs Maximum number of peak pairs to generate
     * @return Vector of peak pairs (anchor peak, target peak)
     */
    static std::vector<std::pair<Peak, Peak>> generate_constellation(
        const std::vector<Peak>& peaks,
        int max_pairs = 0
    );
};

} // namespace audiofp